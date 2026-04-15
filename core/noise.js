// core/noise.js
// CFD v14d: 微小ノイズ項

export function calc_noise(sigma, dt) {
    return sigma * (Math.random() - 0.5) * Math.sqrt(dt);
}