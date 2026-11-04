// Run the app's R code under webR (R compiled to WebAssembly) in Node, as a
// stand-in for the browser: load the data files, run the statistics on the
// fixture list, draw a bullseye to PNG, and print numbers to compare with
// native R (min p and sum of BH per panel).
//
// Setup (once, in a scratch directory):  npm init -y && npm install webr
// Run:  node /path/to/dev/webr_validate.mjs /path/to/app /path/to/output_dir
// Needs network access the first time (webR runtime and the png/jsonlite wasm
// binaries from https://repo.r-wasm.org/).
import { WebR } from 'webr';
import { readFileSync, readdirSync, writeFileSync, statSync } from 'node:fs';
import { join } from 'node:path';

const appDir = process.argv[2];
const webR = new WebR({ repoUrl: 'https://repo.r-wasm.org/' });
await webR.init();
console.log('webR ready');
await webR.installPackages(['png', 'jsonlite'], { quiet: true });
console.log('packages installed');

async function copyDir(src, dst) {
  await webR.FS.mkdir(dst);
  for (const f of readdirSync(src)) {
    const p = join(src, f);
    if (statSync(p).isDirectory()) await copyDir(p, dst + '/' + f);
    else await webR.FS.writeFile(dst + '/' + f, new Uint8Array(readFileSync(p)));
  }
}
await copyDir(appDir, '/app');
console.log('app files copied');

const out = await webR.evalRString(`
  setwd("/app")
  for (f in list.files("R", pattern = "\\\\.R$", full.names = TRUE)) source(f)
  genes <- c("SHANK3","SCN2A","CHD8","SYNGAP1","GRIN2B","MECP2","FMR1","PTEN","TSC1","NRXN1","ADNP","ARID1B","DYRK1A","POGZ","TBR1","ANK2","CNTNAP2","GABRB3","KMT2A","SETD5")
  res <- list()
  for (p in PANELS) {
    t0 <- Sys.time()
    dat <- readRDS(file.path("data", p$file))
    t_load <- as.numeric(Sys.time() - t0, units = "secs")
    t0 <- Sys.time()
    r <- sea_fisher(dat$psi, genes, total = dat$universe$human)
    t_stat <- as.numeric(Sys.time() - t0, units = "secs")
    t05 <- r$table[r$table$threshold == 0.05, ]
    lay <- layout_nodes(p$layout, dat$columns, setNames(t05$n_specific, t05$code))
    t0 <- Sys.time()
    png(paste0("/app/", p$id, ".png"), width = 1000, height = 800)
    bullseye_plot(r$table, lay, tree_dir = ".")
    dev.off()
    t_plot <- as.numeric(Sys.time() - t0, units = "secs")
    res[[p$id]] <- list(found = length(r$found), missing = r$missing, total = r$total,
                        min_p = min(r$table$p), sum_bh = sum(r$table$bh), t_load = t_load, t_stat = t_stat, t_plot = t_plot,
                        png_bytes = file.size(paste0("/app/", p$id, ".png")))
  }
  jsonlite::toJSON(list(r_version = R.version.string, cairo = capabilities("cairo"), panels = res), auto_unbox = TRUE, digits = 12)
`);
console.log(out);
for (const id of ['csea', 'adult', 'dev', 'tsea']) {
  const bytes = await webR.FS.readFile('/app/' + id + '.png');
  writeFileSync(join(process.argv[3], 'webr_' + id + '.png'), bytes);
}
console.log('PNGs written');
process.exit(0);
