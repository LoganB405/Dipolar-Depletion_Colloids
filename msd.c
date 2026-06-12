/*
 * msd.c  —  Compute two MSD datasets from LAMMPS log-spaced dump files.
 *
 * Usage:
 *   ./msd <name> <origin_spacing> <n_origins> <dump_dir> <out_log> <out_linear>
 *
 * Example:
 *   ./msd eps1 10000 10 configs msd_log.dat msd_linear.dat
 *
 * Output format (each file):
 *   # tau   msd   n_samples
 *   ...
 *
 * Compile:
 *   gcc -O2 -o msd msd.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ------------------------------------------------------------------ */
/*  Data structures                                                     */
/* ------------------------------------------------------------------ */

typedef struct {
    double x, y, z;
} Vec3;

typedef struct {
    int   n_atoms;
    Vec3 *pos;
} Snapshot;

typedef struct {
    int       n_frames;
    long     *lag;      /* lag[i] = sim_step[i] - sim_step[0], always >= 0 */
    Snapshot *frames;
} DumpFile;

/* ------------------------------------------------------------------ */
/*  Parsing helpers                                                     */
/* ------------------------------------------------------------------ */

static void skip_line(FILE *f)
{
    int c;
    while ((c = fgetc(f)) != '\n' && c != EOF) {}
}

static int read_frame(FILE *f, Snapshot *snap, long *sim_step)
{
    char line[256];

    /* ITEM: TIMESTEP */
    if (!fgets(line, sizeof(line), f)) return 0;
    if (strncmp(line, "ITEM: TIMESTEP", 14) != 0) return 0;
    if (!fgets(line, sizeof(line), f)) return 0;
    *sim_step = atol(line);

    /* ITEM: NUMBER OF ATOMS */
    if (!fgets(line, sizeof(line), f)) return 0;
    if (!fgets(line, sizeof(line), f)) return 0;
    int n = atoi(line);

    /* ITEM: BOX BOUNDS (3 lines) */
    if (!fgets(line, sizeof(line), f)) return 0;
    skip_line(f); skip_line(f); skip_line(f);

    /* ITEM: ATOMS id type xu yu zu */
    if (!fgets(line, sizeof(line), f)) return 0;

    if (snap->n_atoms != n) {
        snap->pos     = realloc(snap->pos, n * sizeof(Vec3));
        snap->n_atoms = n;
    }

    int    *ids = malloc(n * sizeof(int));
    double *xs  = malloc(n * sizeof(double));
    double *ys  = malloc(n * sizeof(double));
    double *zs  = malloc(n * sizeof(double));

    for (int i = 0; i < n; i++) {
        int    id, type;
        double xu, yu, zu;
        if (!fgets(line, sizeof(line), f)) {
            free(ids); free(xs); free(ys); free(zs);
            return 0;
        }
        sscanf(line, "%d %d %lf %lf %lf", &id, &type, &xu, &yu, &zu);
        ids[i] = id;
        xs[i]  = xu;
        ys[i]  = yu;
        zs[i]  = zu;
    }

    for (int i = 0; i < n; i++) {
        int idx = ids[i] - 1;
        snap->pos[idx].x = xs[i];
        snap->pos[idx].y = ys[i];
        snap->pos[idx].z = zs[i];
    }

    free(ids); free(xs); free(ys); free(zs);
    return 1;
}

/* ------------------------------------------------------------------ */
/*  Load an entire dump file                                            */
/*                                                                      */
/*  Lags are computed as (sim_step - sim_step_of_frame_0), so they     */
/*  are always relative to the first frame in the file regardless of   */
/*  what the global simulation timestep is.                             */
/* ------------------------------------------------------------------ */

static DumpFile *load_dump(const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "Cannot open %s\n", path);
        return NULL;
    }

    DumpFile *df  = calloc(1, sizeof(DumpFile));
    int       cap = 32;
    df->frames    = calloc(cap, sizeof(Snapshot));
    df->lag       = calloc(cap, sizeof(long));

    /* store raw sim_steps first, convert to lags after */
    long *sim_steps = calloc(cap, sizeof(long));

    long sim_step;
    int  frame_idx = 0;

    while (1) {
        if (frame_idx == cap) {
            cap *= 2;
            df->frames  = realloc(df->frames,  cap * sizeof(Snapshot));
            df->lag     = realloc(df->lag,      cap * sizeof(long));
            sim_steps   = realloc(sim_steps,    cap * sizeof(long));
            memset(&df->frames[frame_idx], 0,
                   (cap - frame_idx) * sizeof(Snapshot));
        }

        if (!read_frame(f, &df->frames[frame_idx], &sim_step)) break;

        sim_steps[frame_idx] = sim_step;
        frame_idx++;
    }

    df->n_frames = frame_idx;

    /* Convert absolute sim steps to lags relative to frame 0 */
    long t0 = (frame_idx > 0) ? sim_steps[0] : 0;
    for (int i = 0; i < frame_idx; i++)
        df->lag[i] = sim_steps[i] - t0;

    free(sim_steps);
    fclose(f);
    return df;
}

static void free_dump(DumpFile *df)
{
    for (int i = 0; i < df->n_frames; i++)
        free(df->frames[i].pos);
    free(df->frames);
    free(df->lag);
    free(df);
}

/* ------------------------------------------------------------------ */
/*  MSD between two snapshots                                           */
/* ------------------------------------------------------------------ */

static double msd_between(const Snapshot *a, const Snapshot *b)
{
    double sum = 0.0;
    int n = a->n_atoms;
    for (int i = 0; i < n; i++) {
        double dx = b->pos[i].x - a->pos[i].x;
        double dy = b->pos[i].y - a->pos[i].y;
        double dz = b->pos[i].z - a->pos[i].z;
        sum += dx*dx + dy*dy + dz*dz;
    }
    return sum / n;
}

/* ------------------------------------------------------------------ */
/*  Main                                                                */
/* ------------------------------------------------------------------ */

int main(int argc, char *argv[])
{
    if (argc != 7) {
        fprintf(stderr,
            "Usage: %s <name> <origin_spacing> <n_origins> <dump_dir>"
            " <out_log> <out_linear>\n", argv[0]);
        return 1;
    }

    const char *name           = argv[1];
    long        origin_spacing = atol(argv[2]);
    int         n_origins      = atoi(argv[3]);
    const char *dump_dir       = argv[4];
    const char *out_log_path   = argv[5];
    const char *out_lin_path   = argv[6];

    /* ---- Load all dump files ---- */
    DumpFile **dumps = malloc(n_origins * sizeof(DumpFile *));

    for (int i = 0; i < n_origins; i++) {
        long t0 = (long)i * origin_spacing;
        char path[512];
        snprintf(path, sizeof(path), "%s/dump.%s.%ld.lammpstrj",
                 dump_dir, name, t0);

        printf("Loading %s ...\n", path);
        dumps[i] = load_dump(path);
        if (!dumps[i] || dumps[i]->n_frames == 0) {
            fprintf(stderr, "Failed to load or empty: %s\n", path);
            return 1;
        }
        printf("  -> %d frames, lags:", dumps[i]->n_frames);
        for (int j = 0; j < dumps[i]->n_frames; j++)
            printf(" %ld", dumps[i]->lag[j]);
        printf("\n");
    }

    /* ================================================================
     * DATASET 1: Log-spaced MSD
     *
     * lag[j] is now always relative to frame 0 of that dump file,
     * so tau=0,1,2,4,8,... regardless of global simulation time.
     * ================================================================ */

    FILE *out_log = fopen(out_log_path, "w");
    if (!out_log) { fprintf(stderr, "Cannot open %s\n", out_log_path); return 1; }

    fprintf(out_log, "# LOG-SPACED MSD\n");
    fprintf(out_log, "# tau   msd   n_samples\n");

    /* tau=0: MSD=0 by definition */
    fprintf(out_log, "0  0.000000  %d\n", n_origins);

    int n_log_frames = dumps[0]->n_frames;

    for (int j = 1; j < n_log_frames; j++) {
        long   tau     = dumps[0]->lag[j];
        double msd_sum = 0.0;
        int    count   = 0;

        for (int i = 0; i < n_origins; i++) {
            if (j >= dumps[i]->n_frames) continue;

            if (dumps[i]->lag[j] != tau) {
                fprintf(stderr,
                    "Warning: lag mismatch at origin %d, frame %d: "
                    "expected %ld, got %ld\n",
                    i, j, tau, dumps[i]->lag[j]);
            }

            msd_sum += msd_between(&dumps[i]->frames[0],
                                   &dumps[i]->frames[j]);
            count++;
        }

        if (count > 0)
            fprintf(out_log, "%ld  %.6f  %d\n", tau, msd_sum / count, count);
    }

    fclose(out_log);
    printf("Log-spaced MSD written to %s\n", out_log_path);

    /* ================================================================
     * DATASET 2: Linear-spaced MSD
     *
     * Uses frame[0] from each dump file as the configuration at
     * t = i * origin_spacing.
     * For lag = n * origin_spacing, average over all pairs (i, i+n).
     * ================================================================ */

    FILE *out_lin = fopen(out_lin_path, "w");
    if (!out_lin) { fprintf(stderr, "Cannot open %s\n", out_lin_path); return 1; }

    fprintf(out_lin, "# LINEAR-SPACED MSD\n");
    fprintf(out_lin, "# tau   msd   n_samples\n");

    /* tau=0: MSD=0 by definition */
    fprintf(out_lin, "0  0.000000  %d\n", n_origins);

    for (int n = 1; n < n_origins; n++) {
        long   tau     = (long)n * origin_spacing;
        double msd_sum = 0.0;
        int    count   = 0;

        for (int i = 0; i + n < n_origins; i++) {
            msd_sum += msd_between(&dumps[i]->frames[0],
                                   &dumps[i + n]->frames[0]);
            count++;
        }

        if (count > 0)
            fprintf(out_lin, "%ld  %.6f  %d\n", tau, msd_sum / count, count);
    }

    fclose(out_lin);
    printf("Linear-spaced MSD written to %s\n", out_lin_path);

    /* Cleanup */
    for (int i = 0; i < n_origins; i++) free_dump(dumps[i]);
    free(dumps);

    return 0;
}