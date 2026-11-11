// Text Cube geometry engine (runs in a Web Worker, also testable in node)
// Input masks: Float32Array alpha (R*R, row 0 = top) for front, side, top letters.
const D = 1.1; // half-extent of sampling domain; cube letters live in [-1,1]

function boxBlur(src, R, r) {
  const tmp = new Float32Array(R * R), out = new Float32Array(R * R);
  const w = 2 * r + 1;
  for (let y = 0; y < R; y++) {            // horizontal
    let acc = 0; const row = y * R;
    for (let x = -r; x <= r; x++) acc += src[row + Math.min(R - 1, Math.max(0, x))];
    for (let x = 0; x < R; x++) {
      tmp[row + x] = acc / w;
      acc += src[row + Math.min(R - 1, x + r + 1)] - src[row + Math.max(0, x - r)];
    }
  }
  for (let x = 0; x < R; x++) {            // vertical
    let acc = 0;
    for (let y = -r; y <= r; y++) acc += tmp[Math.min(R - 1, Math.max(0, y)) * R + x];
    for (let y = 0; y < R; y++) {
      out[y * R + x] = acc / w;
      acc += tmp[Math.min(R - 1, y + r + 1) * R + x] - tmp[Math.max(0, y - r) * R + x];
    }
  }
  return out;
}

// Felzenszwalb squared EDT; feature = pixels where seed is true
function edt(seed, R) {
  const INF = 1e20, f = new Float64Array(R), d = new Float64Array(R),
    v = new Int32Array(R), z = new Float64Array(R + 1);
  const g = new Float64Array(R * R);
  for (let i = 0; i < R * R; i++) g[i] = seed[i] ? 0 : INF;
  const pass = (get, set) => {
    for (let q = 0; q < R; q++) f[q] = get(q);
    let k = 0; v[0] = 0; z[0] = -INF; z[1] = INF;
    for (let q = 1; q < R; q++) {
      let s;
      while (true) {
        s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k]);
        if (s <= z[k]) { k--; if (k < 0) { k = 0; break; } } else break;
      }
      k++; v[k] = q; z[k] = s; z[k + 1] = INF;
    }
    k = 0;
    for (let q = 0; q < R; q++) {
      while (z[k + 1] < q) k++;
      d[q] = (q - v[k]) * (q - v[k]) + f[v[k]];
    }
    for (let q = 0; q < R; q++) set(q, d[q]);
  };
  for (let x = 0; x < R; x++) pass(q => g[q * R + x], (q, val) => { g[q * R + x] = val; });
  for (let y = 0; y < R; y++) pass(q => g[y * R + q], (q, val) => { g[y * R + q] = val; });
  return g;
}

function makeSDF(alpha, R, roundUnits, boldUnits) {
  const px = 2 * D / R;
  let a = alpha;
  const rpx = roundUnits / px;
  if (rpx >= 0.75) {
    // three box passes ≈ gaussian; radius chosen so corner radius ≈ roundUnits
    const r = Math.max(1, Math.round(rpx * 0.55));
    a = boxBlur(boxBlur(boxBlur(a, R, r), R, r), R, r);
  }
  const inside = new Uint8Array(R * R), outside = new Uint8Array(R * R);
  for (let i = 0; i < R * R; i++) { inside[i] = a[i] >= 0.5 ? 1 : 0; outside[i] = 1 - inside[i]; }
  const dToIn = edt(inside, R), dToOut = edt(outside, R);
  const s = new Float32Array(R * R);
  for (let i = 0; i < R * R; i++) {
    const d = inside[i] ? -(Math.sqrt(dToOut[i]) - 0.5) : (Math.sqrt(dToIn[i]) - 0.5);
    s[i] = d * px - boldUnits;
  }
  return s;
}

// bilinear sample; u right, v up, domain coords
function sample(s, R, u, v) {
  let fx = (u + D) / (2 * D) * R - 0.5, fy = (D - v) / (2 * D) * R - 0.5;
  let extra = 0;
  const cx = Math.min(R - 1, Math.max(0, fx)), cy = Math.min(R - 1, Math.max(0, fy));
  extra = Math.hypot(fx - cx, fy - cy) * (2 * D / R);
  const x0 = Math.min(R - 2, Math.floor(cx)), y0 = Math.min(R - 2, Math.floor(cy));
  const tx = cx - x0, ty = cy - y0, i = y0 * R + x0;
  const a = s[i] * (1 - tx) + s[i + 1] * tx, b = s[i + R] * (1 - tx) + s[i + R + 1] * tx;
  return a * (1 - ty) + b * ty + extra;
}

function smax(a, b, k) {
  if (k <= 0) return a > b ? a : b;
  const h = Math.max(k - Math.abs(a - b), 0) / k;
  return (a > b ? a : b) + h * h * k * 0.25;
}

function build(p) {
  const t0 = Date.now();
  const R = p.R, N = p.N, h = 2 * D / (N - 1);
  const pos = i => -D + i * h;
  const sF = makeSDF(p.masks[0], R, p.round, p.bold);
  const sE = makeSDF(p.masks[1], R, p.round, p.bold);
  const sW = makeSDF(p.masks[2], R, p.round, p.bold);
  // tables: TF[i + N*j] = F(x_i, y_j); TE[k + N*j] = E(-z_k, y_j); TW[i + N*k] = W(x_i, -z_k)
  const TF = new Float32Array(N * N), TE = new Float32Array(N * N), TW = new Float32Array(N * N);
  for (let a = 0; a < N; a++) for (let b = 0; b < N; b++) {
    TF[a + N * b] = sample(sF, R, pos(a), pos(b));
    TE[a + N * b] = sample(sE, R, -pos(a), pos(b));
    TW[a + N * b] = sample(sW, R, pos(a), -pos(b));
  }
  const NN = N * N, f = new Float32Array(N * NN);
  const k = p.fillet;
  for (let kk = 0; kk < N; kk++) for (let j = 0; j < N; j++) {
    const e = TE[kk + N * j], base = N * (j + N * kk);
    for (let i = 0; i < N; i++) {
      f[base + i] = smax(smax(TF[i + N * j], e, k), TW[i + N * kk], k);
    }
  }
  // force boundary outside so the mesh is closed
  for (let kk = 0; kk < N; kk++) for (let j = 0; j < N; j++) for (let i = 0; i < N; i++)
    if (i === 0 || j === 0 || kk === 0 || i === N - 1 || j === N - 1 || kk === N - 1) {
      const id = i + N * (j + N * kk); if (f[id] < h) f[id] = h;
    }

  // connected components of solid voxels (6-connectivity)
  let pieces = 0, removed = 0;
  {
    const lab = new Int32Array(N * NN), q = new Int32Array(N * NN);
    const sizes = [0];
    for (let s = 0; s < N * NN; s++) {
      if (f[s] >= 0 || lab[s]) continue;
      const L = sizes.length; let head = 0, tail = 0, cnt = 0;
      q[tail++] = s; lab[s] = L;
      while (head < tail) {
        const c = q[head++]; cnt++;
        const i = c % N, j = ((c / N) | 0) % N, kk = (c / NN) | 0;
        if (i > 0) { const n = c - 1; if (!lab[n] && f[n] < 0) { lab[n] = L; q[tail++] = n; } }
        if (i < N - 1) { const n = c + 1; if (!lab[n] && f[n] < 0) { lab[n] = L; q[tail++] = n; } }
        if (j > 0) { const n = c - N; if (!lab[n] && f[n] < 0) { lab[n] = L; q[tail++] = n; } }
        if (j < N - 1) { const n = c + N; if (!lab[n] && f[n] < 0) { lab[n] = L; q[tail++] = n; } }
        if (kk > 0) { const n = c - NN; if (!lab[n] && f[n] < 0) { lab[n] = L; q[tail++] = n; } }
        if (kk < N - 1) { const n = c + NN; if (!lab[n] && f[n] < 0) { lab[n] = L; q[tail++] = n; } }
      }
      sizes.push(cnt);
    }
    let best = 0, bestL = 0;
    for (let L = 1; L < sizes.length; L++) { if (sizes[L] > best) { best = sizes[L]; bestL = L; } }
    pieces = sizes.length - 1;
    if (p.keepLargest && pieces > 1) {
      for (let s = 0; s < N * NN; s++) if (lab[s] && lab[s] !== bestL) { f[s] = h * 0.5; removed++; }
    }
  }

  // silhouette coverage per view. code: 0 empty, 1 ok, 2 missing, 3 extra
  const views = [];
  const projF = new Uint8Array(NN), projE = new Uint8Array(NN), projW = new Uint8Array(NN);
  for (let kk = 0; kk < N; kk++) for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
    if (f[i + N * (j + N * kk)] < 0) { projF[i + N * j] = 1; projE[kk + N * j] = 1; projW[i + N * kk] = 1; }
  }
  // image layout for thumbnails: row 0 = top of view, col 0 = left of view
  const mk = (T, P, toImg) => {
    const img = new Uint8Array(NN); let tgt = 0, ok = 0;
    for (let b = 0; b < N; b++) for (let a = 0; a < N; a++) {
      const want = T[a + N * b] < -h, has = P[a + N * b], loose = T[a + N * b] < h;
      let c = 0;
      if (want) { tgt++; if (has) { ok++; c = 1; } else c = 2; }
      else if (has) c = loose ? 1 : 3;
      img[toImg(a, b)] = c;
    }
    return { img, coverage: tgt ? ok / tgt : 1 };
  };
  // front: a=i (x right), b=j (y up) → row = N-1-j
  views.push(mk(TF, projF, (a, b) => (N - 1 - b) * N + a));
  // side: a=k, screen right = -z → col = N-1-k, row = N-1-j
  views.push(mk(TE, projE, (a, b) => (N - 1 - b) * N + (N - 1 - a)));
  // top: a=i (x right), b=k; screen up = -z → top row is smallest z → row = k
  views.push(mk(TW, projW, (a, b) => b * N + a));

  // ---- surface nets ----
  const M = N - 1;
  const cellV = new Int32Array(M * M * M).fill(-1);
  const P = [];
  const corner = [[0,0,0],[1,0,0],[0,1,0],[1,1,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1]];
  const edges = [[0,1],[2,3],[4,5],[6,7],[0,2],[1,3],[4,6],[5,7],[0,4],[1,5],[2,6],[3,7]];
  const val = new Float32Array(8);
  for (let kk = 0; kk < M; kk++) for (let j = 0; j < M; j++) for (let i = 0; i < M; i++) {
    let mask = 0;
    for (let c = 0; c < 8; c++) {
      const o = corner[c]; const v = f[(i + o[0]) + N * ((j + o[1]) + N * (kk + o[2]))];
      val[c] = v; if (v < 0) mask |= 1 << c;
    }
    if (mask === 0 || mask === 255) continue;
    let sx = 0, sy = 0, sz = 0, n = 0;
    for (const [a, b] of edges) {
      const va = val[a], vb = val[b];
      if ((va < 0) === (vb < 0)) continue;
      const t = va / (va - vb), A = corner[a], B = corner[b];
      sx += A[0] + (B[0] - A[0]) * t; sy += A[1] + (B[1] - A[1]) * t; sz += A[2] + (B[2] - A[2]) * t; n++;
    }
    cellV[i + M * (j + M * kk)] = P.length / 3;
    P.push(pos(i) + h * sx / n, pos(j) + h * sy / n, pos(kk) + h * sz / n);
  }
  const Q = []; // quads
  const cv = (i, j, kk) => cellV[i + M * (j + M * kk)];
  const quad = (a, b, c, d, flip) => { if (flip) Q.push(d, c, b, a); else Q.push(a, b, c, d); };
  for (let kk = 1; kk < M; kk++) for (let j = 1; j < M; j++) for (let i = 0; i < M; i++) {
    const v0 = f[i + N * (j + N * kk)], v1 = f[i + 1 + N * (j + N * kk)];
    if ((v0 < 0) !== (v1 < 0)) quad(cv(i, j-1, kk-1), cv(i, j, kk-1), cv(i, j, kk), cv(i, j-1, kk), !(v0 < 0));
  }
  for (let kk = 1; kk < M; kk++) for (let j = 0; j < M; j++) for (let i = 1; i < M; i++) {
    const v0 = f[i + N * (j + N * kk)], v1 = f[i + N * (j + 1 + N * kk)];
    if ((v0 < 0) !== (v1 < 0)) quad(cv(i-1, j, kk-1), cv(i-1, j, kk), cv(i, j, kk), cv(i, j, kk-1), !(v0 < 0));
  }
  for (let kk = 0; kk < M; kk++) for (let j = 1; j < M; j++) for (let i = 1; i < M; i++) {
    const v0 = f[i + N * (j + N * kk)], v1 = f[i + N * (j + N * (kk + 1))];
    if ((v0 < 0) !== (v1 < 0)) quad(cv(i-1, j-1, kk), cv(i, j-1, kk), cv(i, j, kk), cv(i-1, j, kk), !(v0 < 0));
  }
  const nv = P.length / 3, positions = new Float32Array(P);

  // ---- Taubin smoothing (polish) ----
  if (p.polish > 0 && nv) {
    const deg = new Int32Array(nv + 1);
    for (let q = 0; q < Q.length; q += 4) for (let e = 0; e < 4; e++) { deg[Q[q + e]]++; deg[Q[q + ((e + 1) & 3)]]++; }
    const start = new Int32Array(nv + 1);
    for (let v = 0; v < nv; v++) start[v + 1] = start[v] + deg[v];
    const fillp = start.slice(0, nv), adj = new Int32Array(start[nv]);
    for (let q = 0; q < Q.length; q += 4) for (let e = 0; e < 4; e++) {
      const a = Q[q + e], b = Q[q + ((e + 1) & 3)];
      adj[fillp[a]++] = b; adj[fillp[b]++] = a;
    }
    const tmp = new Float32Array(nv * 3);
    const step = lam => {
      for (let v = 0; v < nv; v++) {
        let ax = 0, ay = 0, az = 0; const s = start[v], e = start[v + 1], c = e - s;
        if (!c) { tmp[3*v] = positions[3*v]; tmp[3*v+1] = positions[3*v+1]; tmp[3*v+2] = positions[3*v+2]; continue; }
        for (let t = s; t < e; t++) { const u = adj[t] * 3; ax += positions[u]; ay += positions[u + 1]; az += positions[u + 2]; }
        tmp[3*v] = positions[3*v] + lam * (ax / c - positions[3*v]);
        tmp[3*v+1] = positions[3*v+1] + lam * (ay / c - positions[3*v+1]);
        tmp[3*v+2] = positions[3*v+2] + lam * (az / c - positions[3*v+2]);
      }
      positions.set(tmp);
    };
    for (let it = 0; it < p.polish; it++) { step(0.5); step(-0.53); }
  }

  // triangulate quads along the shorter diagonal
  const idx = new Uint32Array(Q.length / 4 * 6);
  const d2 = (a, b) => { const x = positions[3*a] - positions[3*b], y = positions[3*a+1] - positions[3*b+1], z = positions[3*a+2] - positions[3*b+2]; return x*x + y*y + z*z; };
  for (let q = 0, t = 0; q < Q.length; q += 4) {
    const a = Q[q], b = Q[q + 1], c = Q[q + 2], d = Q[q + 3];
    if (d2(a, c) <= d2(b, d)) { idx[t++] = a; idx[t++] = b; idx[t++] = c; idx[t++] = a; idx[t++] = c; idx[t++] = d; }
    else { idx[t++] = a; idx[t++] = b; idx[t++] = d; idx[t++] = b; idx[t++] = c; idx[t++] = d; }
  }
  return { positions, index: idx, N, views, pieces, removed, ms: Date.now() - t0, id: p.id };
}

if (typeof self !== 'undefined' && typeof window === 'undefined' && typeof module === 'undefined') {
  self.onmessage = e => {
    const r = build(e.data);
    self.postMessage(r, [r.positions.buffer, r.index.buffer, ...r.views.map(v => v.img.buffer)]);
  };
}
if (typeof module !== 'undefined') module.exports = { build };
