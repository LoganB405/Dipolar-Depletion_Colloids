/*
 * gr.c  –  Radial distribution function g(r)
 *
 * Averages over all frames in all LAMMPS trajectory files (.lammpstrj)
 * found in a given directory.
 *
 * Expected dump format:
 *   dump d all custom ... id type xu yu zu
 *
 * Unwrapped coordinates are wrapped back into the primary image before
 * computing pair distances (minimum-image PBC).
 *
 * Usage:
 *   gcc -O2 gr.c -lm -o gr
 *   ./gr <trajectory_dir> <dr> <r_max> [output_file]
 *
 *   trajectory_dir  – folder containing *.lammpstrj files
 *   dr              – bin width (LJ units, e.g. 0.02)
 *   r_max           – maximum r; should be <= L/2
 *   output_file     – optional; defaults to "gr_output.dat"
 *
 * Output columns:  r   g(r)
 *
 * Normalization (exact for fixed N, V):
 *
 *   g(r) = hist[b] / (total_frames * N * rho * V_shell)
 *
 * where:
 *   hist[b]       counts BOTH i→j and j→i (factor of 2 already in)
 *   N             true atom count per frame (tracked explicitly)
 *   rho           = N / V
 *   V_shell       = (4/3)*pi*(r_hi^3 - r_lo^3)
 *
 * Equivalently, dividing by N cancels and you get:
 *   ideal = total_frames * (N-1) * rho * V_shell
 * which is used below.
 */

 #define _GNU_SOURCE
 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <math.h>
 #include <dirent.h>
 #include <errno.h>
 
 /* ------------------------------------------------------------------ */
 /* Data structures                                                      */
 /* ------------------------------------------------------------------ */
 
 typedef struct { double x, y, z; } Vec3;
 
 typedef struct {
     double xlo, xhi;
     double ylo, yhi;
     double zlo, zhi;
 } Box;
 
 /* ------------------------------------------------------------------ */
 /* Helpers                                                              */
 /* ------------------------------------------------------------------ */
 
 static inline double pbc_wrap(double dx, double L)
 {
     dx -= L * round(dx / L);
     return dx;
 }
 
 static inline double dist2_pbc(const Vec3 *a, const Vec3 *b, const Box *box)
 {
     double dx = pbc_wrap(a->x - b->x, box->xhi - box->xlo);
     double dy = pbc_wrap(a->y - b->y, box->yhi - box->ylo);
     double dz = pbc_wrap(a->z - b->z, box->zhi - box->zlo);
     return dx*dx + dy*dy + dz*dz;
 }
 
 /* ------------------------------------------------------------------ */
 /* Frame reader                                                          */
 /* ------------------------------------------------------------------ */
 
 static int read_frame(FILE *fp, Vec3 **atoms, int *n_atoms, Box *box)
 {
     char line[512];
     int  n = 0;
 
     if (!fgets(line, sizeof(line), fp)) return 0;
     if (strncmp(line, "ITEM: TIMESTEP", 14) != 0) return 0;
     if (!fgets(line, sizeof(line), fp)) return 0;
 
     if (!fgets(line, sizeof(line), fp)) return 0;
     if (strncmp(line, "ITEM: NUMBER OF ATOMS", 21) != 0) return 0;
     if (!fgets(line, sizeof(line), fp)) return 0;
     n = atoi(line);
     if (n <= 0) return 0;
 
     if (!fgets(line, sizeof(line), fp)) return 0;
     if (strncmp(line, "ITEM: BOX BOUNDS", 16) != 0) return 0;
     if (!fgets(line, sizeof(line), fp)) return 0;
     sscanf(line, "%lf %lf", &box->xlo, &box->xhi);
     if (!fgets(line, sizeof(line), fp)) return 0;
     sscanf(line, "%lf %lf", &box->ylo, &box->yhi);
     if (!fgets(line, sizeof(line), fp)) return 0;
     sscanf(line, "%lf %lf", &box->zlo, &box->zhi);
 
     if (!fgets(line, sizeof(line), fp)) return 0;
     if (strncmp(line, "ITEM: ATOMS", 11) != 0) return 0;
 
     /* Detect column indices for xu, yu, zu */
     int col_xu = -1, col_yu = -1, col_zu = -1;
     {
         char header[512];
         strncpy(header, line + 11, sizeof(header) - 1);
         header[sizeof(header)-1] = '\0';
         char *tok = strtok(header, " \t\r\n");
         int col = 0;
         while (tok) {
             if      (strcmp(tok, "xu") == 0) col_xu = col;
             else if (strcmp(tok, "yu") == 0) col_yu = col;
             else if (strcmp(tok, "zu") == 0) col_zu = col;
             tok = strtok(NULL, " \t\r\n");
             col++;
         }
     }
     if (col_xu < 0) col_xu = 2;
     if (col_yu < 0) col_yu = 3;
     if (col_zu < 0) col_zu = 4;
     int max_col = col_xu > col_yu ? col_xu : col_yu;
     if (col_zu > max_col) max_col = col_zu;
 
     Vec3 *arr = realloc(*atoms, n * sizeof(Vec3));
     if (!arr) { fprintf(stderr, "Out of memory.\n"); exit(1); }
     *atoms   = arr;
     *n_atoms = n;
 
     double lx = box->xhi - box->xlo;
     double ly = box->yhi - box->ylo;
     double lz = box->zhi - box->zlo;
 
     for (int i = 0; i < n; i++) {
         if (!fgets(line, sizeof(line), fp)) return 0;
         char buf[512];
         strncpy(buf, line, sizeof(buf) - 1);
         buf[sizeof(buf)-1] = '\0';
 
         double vals[16] = {0};
         char *tok = strtok(buf, " \t\r\n");
         int col = 0;
         while (tok && col <= max_col && col < 16) {
             vals[col] = atof(tok);
             tok = strtok(NULL, " \t\r\n");
             col++;
         }
 
         double xu = vals[col_xu];
         double yu = vals[col_yu];
         double zu = vals[col_zu];
 
         /* Wrap unwrapped coordinates into primary image */
         xu -= lx * floor((xu - box->xlo) / lx);
         yu -= ly * floor((yu - box->ylo) / ly);
         zu -= lz * floor((zu - box->zlo) / lz);
 
         arr[i].x = xu;
         arr[i].y = yu;
         arr[i].z = zu;
     }
 
     return 1;
 }
 
 /* ------------------------------------------------------------------ */
 /* Accumulate histogram from one frame                                  */
 /* ------------------------------------------------------------------ */
 
 static void accumulate_frame(const Vec3 *atoms, int n, const Box *box,
                              double dr, int n_bins, long long *hist)
 {
     double r2_max = (n_bins * dr) * (n_bins * dr);
 
     for (int i = 0; i < n - 1; i++) {
         for (int j = i + 1; j < n; j++) {
             double r2 = dist2_pbc(&atoms[i], &atoms[j], box);
             if (r2 > 0.0 && r2 < r2_max) {
                 int b = (int)(sqrt(r2) / dr);
                 if (b < n_bins)
                     hist[b] += 2;   /* count i→j and j→i */
             }
         }
     }
 }
 
 /* ------------------------------------------------------------------ */
 /* Process one trajectory file                                          */
 /* ------------------------------------------------------------------ */
 
 static void process_file(const char *path,
                          double dr, int n_bins,
                          long long *hist,
                          long long *total_frames,
                          double    *sum_rho,
                          long long *sum_N)        /* <-- explicit N tracking */
 {
     FILE *fp = fopen(path, "r");
     if (!fp) {
         fprintf(stderr, "Warning: cannot open %s: %s\n", path, strerror(errno));
         return;
     }
 
     Vec3 *atoms   = NULL;
     int   n_atoms = 0;
     Box   box;
     long long frames_this_file = 0;
 
     while (read_frame(fp, &atoms, &n_atoms, &box)) {
         double vol = (box.xhi - box.xlo) *
                      (box.yhi - box.ylo) *
                      (box.zhi - box.zlo);
 
         *sum_rho += (double)n_atoms / vol;
         *sum_N   += n_atoms;                /* track N directly */
 
         accumulate_frame(atoms, n_atoms, &box, dr, n_bins, hist);
         frames_this_file++;
     }
 
     *total_frames += frames_this_file;
     printf("  %s : %lld frames\n", path, frames_this_file);
 
     free(atoms);
     fclose(fp);
 }
 
 /* ------------------------------------------------------------------ */
 /* main                                                                 */
 /* ------------------------------------------------------------------ */
 
 int main(int argc, char *argv[])
 {
     if (argc < 4) {
         fprintf(stderr,
             "Usage: %s <trajectory_dir> <dr> <r_max> [output_file]\n",
             argv[0]);
         return 1;
     }
 
     const char *traj_dir = argv[1];
     double      dr       = atof(argv[2]);
     double      r_max    = atof(argv[3]);
     const char *out_file = (argc >= 5) ? argv[4] : "gr_output.dat";
 
     if (dr <= 0.0 || r_max <= 0.0) {
         fprintf(stderr, "dr and r_max must be positive.\n");
         return 1;
     }
 
     int n_bins = (int)ceil(r_max / dr);
     printf("Bins: %d  (dr=%.4f, r_max=%.4f)\n", n_bins, dr, r_max);
 
     long long *hist = calloc(n_bins, sizeof(long long));
     if (!hist) { fprintf(stderr, "Out of memory.\n"); return 1; }
 
     long long total_frames = 0;
     long long sum_N        = 0;
     double    sum_rho      = 0.0;
     int       n_files      = 0;
 
     /* ---- Iterate over *.lammpstrj files ---- */
     DIR *dir = opendir(traj_dir);
     if (!dir) {
         fprintf(stderr, "Cannot open directory '%s': %s\n",
                 traj_dir, strerror(errno));
         free(hist);
         return 1;
     }
 
     struct dirent *ent;
     while ((ent = readdir(dir)) != NULL) {
         const char *name = ent->d_name;
         size_t len = strlen(name);
         if (len < 10) continue;
         if (strcmp(name + len - 10, ".lammpstrj") != 0) continue;
 
         char path[1024];
         snprintf(path, sizeof(path), "%s/%s", traj_dir, name);
         n_files++;
         process_file(path, dr, n_bins, hist, &total_frames, &sum_rho, &sum_N);
     }
     closedir(dir);
 
     if (n_files == 0) {
         fprintf(stderr, "No .lammpstrj files found in '%s'.\n", traj_dir);
         free(hist);
         return 1;
     }
     if (total_frames == 0) {
         fprintf(stderr, "No frames read.\n");
         free(hist);
         return 1;
     }
 
     /* True averages — no estimation tricks */
     double N_avg  = (double)sum_N  / (double)total_frames;
     double rho    = sum_rho / (double)total_frames;
 
     /* Warn if r_max looks larger than L/2 (approximate check) */
     double V_avg  = N_avg / rho;
     double L_half = 0.5 * cbrt(V_avg);
     if (r_max > L_half) {
         fprintf(stderr,
             "Warning: r_max (%.3f) > L/2 (%.3f). "
             "Minimum-image convention may be violated.\n",
             r_max, L_half);
     }
 
     printf("Total files : %d\n",      n_files);
     printf("Total frames: %lld\n",    total_frames);
     printf("Mean N      : %.1f\n",    N_avg);
     printf("Mean rho    : %.6f\n",    rho);
 
     /* ---- Normalise to g(r) ---- */
     /*
      * hist[b] counts ordered pairs (i→j AND j→i) landing in bin b,
      * summed over all frames.
      *
      * For a uniform ideal gas the expected count is:
      *   ideal = total_frames * N * rho * V_shell
      *
      * Because hist already has the factor of 2 (both directions),
      * this is equivalent to:
      *   ideal = total_frames * (N-1) * rho * V_shell
      *
      * The two forms differ only by N vs (N-1), which matters for
      * small systems. We use (N-1) for correctness.
      *
      * Result: g(r) → 1 at large r for any density.
      */
     FILE *fout = fopen(out_file, "w");
     if (!fout) {
         fprintf(stderr, "Cannot open output '%s': %s\n",
                 out_file, strerror(errno));
         free(hist);
         return 1;
     }
 
     fprintf(fout, "# g(r) averaged over %d files, %lld frames\n",
             n_files, total_frames);
     fprintf(fout, "# Mean N   = %.1f\n",    N_avg);
     fprintf(fout, "# Mean rho = %.6f\n",    rho);
     fprintf(fout, "# r   g(r)\n");
 
     for (int b = 0; b < n_bins; b++) {
         double r_lo  = b * dr;
         double r_hi  = r_lo + dr;
         double r_mid = r_lo + 0.5 * dr;
 
         double v_shell = (4.0 / 3.0) * M_PI *
                          (r_hi*r_hi*r_hi - r_lo*r_lo*r_lo);
 
         /* Expected ideal-gas count for ordered pairs */
         double ideal = (double)total_frames * (N_avg - 1.0) * rho * v_shell;
 
         double gr = (ideal > 0.0) ? (double)hist[b] / ideal : 0.0;
 
         fprintf(fout, "%.6f  %.6f\n", r_mid, gr);
     }
 
     fclose(fout);
     free(hist);
 
     printf("g(r) written to '%s'\n", out_file);
     return 0;
 }
