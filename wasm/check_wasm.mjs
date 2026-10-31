// The browser build against the certified grid, and against native results if given:
//   node check_wasm.mjs <eliah_wasm.wasm> <eckart_certified.csv> [native.json]
import { readFileSync } from 'node:fs';

const [wasmPath, csvPath, nativePath] = process.argv.slice(2);
const { instance } = await WebAssembly.instantiate(readFileSync(wasmPath), {});
const e = instance.exports;
const out = () => new Float64Array(e.memory.buffer, e.result_ptr(), 6);
const lines = readFileSync(csvPath, 'utf8').trim().split(/\r?\n/);
const head = lines[0].split(',');
const ix = (n) => head.indexOf(n);
const native = nativePath ? JSON.parse(readFileSync(nativePath, 'utf8')) : null;
let n = 0, bad = 0, worst = 0, cosWrong = 0, identical = 0, worstNative = 0;
for (const line of lines.slice(1)) {
  const f = line.split(',');
  const v1 = Number(f[ix('v1')]), v2 = Number(f[ix('v2')]), u = Number(f[ix('ustar')]);
  const ref = Number(f[ix('kappa_mid')]), cos = f[ix('cos_branch')] === '1';
  const rc = e.kappa(v1, v2, u);
  const r = out();
  const rel = rc === 0 && Number.isFinite(r[0]) ? Math.abs(r[0] / ref - 1) : Infinity;
  worst = Math.max(worst, rel);
  if (rel > 1e-12) bad++;
  if ((r[4] === 1) !== cos) cosWrong++;
  if (native) {
    const nv = native[n];
    if (Object.is(nv, r[0])) identical++;
    worstNative = Math.max(worstNative, Math.abs(r[0] / nv - 1));
  }
  n++;
}
console.log(`${n} certified points: ${n - bad} within 1e-12, worst relative error ${worst.toExponential(2)}, ` +
            `cos-branch flag mismatches ${cosWrong}`);
if (native) console.log(`against the native library: ${identical} of ${n} bit-identical, worst relative difference ${worstNative.toExponential(2)}`);
const badInput = e.kappa(0, 1, 1);
console.log(`invalid input returns ${badInput} (expected 1)`);
process.exit(n === 1197 && bad === 0 && cosWrong === 0 && badInput === 1 ? 0 : 1);
