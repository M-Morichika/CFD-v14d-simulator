// /viz/render.js
// CFD v14d visualization renderer

export class Renderer {
    constructor(heatCanvasId, chartCanvasId, loopCanvasId, N) {
        this.heatCtx = document.getElementById(heatCanvasId).getContext('2d');
        this.chartCtx = document.getElementById(chartCanvasId).getContext('2d');
        this.loopCtx = document.getElementById(loopCanvasId).getContext('2d');
        this.N = N;
        this.width = this.heatCtx.canvas.width;
        this.height = this.heatCtx.canvas.height;
    }

    drawHeatmap(rho_array, F_array, max_rho, betaTrail) {
        let maxColorVal = Math.max(0.1, max_rho);
        if (isNaN(maxColorVal) || !isFinite(maxColorVal)) maxColorVal = 1.0;

        let imgData = this.heatCtx.createImageData(this.N, this.N);

        for (let i = 0; i < this.N; i++) {
            for (let j = 0; j < this.N; j++) {
                let idx = (i * this.N + j) * 4;
                let val = rho_array[i][j] / maxColorVal;

                // Viridis-like color mapping
                let r = Math.min(255, Math.max(0, 255 * (val - 0.2) * 2));
                let g = Math.min(255, Math.max(0, 255 * val));
                let b = Math.min(255, Math.max(0, 255 * (1 - val * 0.5) + 50));

                imgData.data[idx] = r; imgData.data[idx + 1] = g;
                imgData.data[idx + 2] = b; imgData.data[idx + 3] = 255;
            }
        }

        let tempCanvas = document.createElement('canvas');
        tempCanvas.width = this.N; tempCanvas.height = this.N;
        tempCanvas.getContext('2d').putImageData(imgData, 0, 0);

        this.heatCtx.imageSmoothingEnabled = false;
        this.heatCtx.drawImage(tempCanvas, 0, 0, this.width, this.height);

        // F Contours (ポテンシャル地形の等高線)
        if (F_array) {
            this.heatCtx.strokeStyle = "rgba(255, 255, 255, 0.25)";
            this.heatCtx.lineWidth = 1;
            let contourInterval = 0.05;

            for (let i = 1; i < this.N - 1; i++) {
                for (let j = 1; j < this.N - 1; j++) {
                    let f_val = F_array[i][j];
                    if (Math.floor(f_val / contourInterval) !== Math.floor(F_array[i][j + 1] / contourInterval) ||
                        Math.floor(f_val / contourInterval) !== Math.floor(F_array[i + 1][j] / contourInterval)) {
                        let px = (j / this.N) * this.width;
                        let py = (i / this.N) * this.height;
                        this.heatCtx.strokeRect(px, py, 1, 1);
                    }
                }
            }
        }

        // β* Trajectory (注意の焦点の軌跡)
        if (betaTrail && betaTrail.length > 1) {
            this.heatCtx.beginPath();
            this.heatCtx.strokeStyle = "rgba(255, 200, 0, 0.8)";
            this.heatCtx.lineWidth = 2;
            for (let k = 0; k < betaTrail.length; k++) {
                let px = ((betaTrail[k].j + 0.5) / this.N) * this.width;
                let py = ((betaTrail[k].i + 0.5) / this.N) * this.height;
                if (k === 0) this.heatCtx.moveTo(px, py);
                else this.heatCtx.lineTo(px, py);
            }
            this.heatCtx.stroke();

            let last = betaTrail[betaTrail.length - 1];
            let px = ((last.j + 0.5) / this.N) * this.width;
            let py = ((last.i + 0.5) / this.N) * this.height;
            this.heatCtx.beginPath();
            this.heatCtx.arc(px, py, 4, 0, Math.PI * 2);
            this.heatCtx.fillStyle = "white";
            this.heatCtx.fill();
        }
    }

    drawChart(history_M, history_pi, history_L) {
        let w = this.chartCtx.canvas.width;
        let h = this.chartCtx.canvas.height;
        this.chartCtx.clearRect(0, 0, w, h);

        let h_top = h * 0.6;
        let h_bot = h * 0.4;
        let stepX = w / history_M.length;

        // 区切り線
        this.chartCtx.strokeStyle = "#444";
        this.chartCtx.lineWidth = 1;
        this.chartCtx.beginPath();
        this.chartCtx.moveTo(0, h_top);
        this.chartCtx.lineTo(w, h_top);
        this.chartCtx.stroke();

        let validM = history_M.filter(v => !isNaN(v));
        let maxM = validM.length > 0 ? Math.max(...validM, 0.1) : 0.1;
        let scaleM = maxM * 1.2; // 動的スケーリング
        let scaleL = 1.2;

        // M(t)
        this.chartCtx.beginPath();
        this.chartCtx.strokeStyle = "#4a90e2";
        this.chartCtx.lineWidth = 2;
        for (let i = 0; i < history_M.length; i++) {
            let x = i * stepX;
            let y = h_top - (history_M[i] / scaleM) * h_top;
            if (i === 0) this.chartCtx.moveTo(x, y);
            else this.chartCtx.lineTo(x, y);
        }
        this.chartCtx.stroke();

        // L(t) (フォーカスの局在度合い)
        if (history_L) {
            this.chartCtx.beginPath();
            this.chartCtx.strokeStyle = "rgba(255, 200, 0, 1.0)";
            this.chartCtx.lineWidth = 2;
            for (let i = 0; i < history_L.length; i++) {
                let x = i * stepX;
                let y = h_top - (history_L[i] / scaleL) * h_top;
                if (i === 0) this.chartCtx.moveTo(x, y);
                else this.chartCtx.lineTo(x, y);
            }
            this.chartCtx.stroke();
        }

        // pi_int(t)
        if (history_pi) {
            let scalePi = 3.0; // v14d の動的レンジ (0.0 ~ 3.0) に合わせる
            this.chartCtx.beginPath();
            this.chartCtx.strokeStyle = "#ff5555";
            this.chartCtx.lineWidth = 2;
            for (let i = 0; i < history_pi.length; i++) {
                let x = i * stepX;
                let y = h_top + h_bot - (history_pi[i] / scalePi) * h_bot;
                if (i === 0) this.chartCtx.moveTo(x, y);
                else this.chartCtx.lineTo(x, y);
            }
            this.chartCtx.stroke();

            // λ₀ (1.0) の点火閾値ライン
            let y_thresh = h_top + h_bot - (1.0 / scalePi) * h_bot;
            this.chartCtx.strokeStyle = "#888";
            this.chartCtx.setLineDash([3, 3]);
            this.chartCtx.beginPath();
            this.chartCtx.moveTo(0, y_thresh);
            this.chartCtx.lineTo(w, y_thresh);
            this.chartCtx.stroke();
            this.chartCtx.setLineDash([]);
        }

        this.chartCtx.font = "10px sans-serif";
        this.chartCtx.fillStyle = "#4a90e2";
        this.chartCtx.fillText(`M(t)`, 5, 12);
        this.chartCtx.fillStyle = "#ff5555";
        this.chartCtx.fillText(`π_int(t)`, 5, h_top + 12);
    }

    drawLoop(history_pi, history_M) {
        let w = this.loopCtx.canvas.width;
        let h = this.loopCtx.canvas.height;
        this.loopCtx.clearRect(0, 0, w, h);

        let padX = 25;
        let padY = 20;

        // 軸の描画
        this.loopCtx.strokeStyle = "#444";
        this.loopCtx.lineWidth = 1;
        this.loopCtx.beginPath();
        this.loopCtx.moveTo(padX, h - padY);
        this.loopCtx.lineTo(w - 10, h - padY);
        this.loopCtx.moveTo(padX, 10);
        this.loopCtx.lineTo(padX, h - padY);
        this.loopCtx.stroke();

        this.loopCtx.font = "10px sans-serif";
        this.loopCtx.fillStyle = "#aaa";
        this.loopCtx.fillText("π_int", w / 2, h - 5);
        this.loopCtx.save();
        this.loopCtx.translate(10, h / 2 + 10);
        this.loopCtx.rotate(-Math.PI / 2);
        this.loopCtx.fillText("M", 0, 0);
        this.loopCtx.restore();

        // v14d向け: Mの動的スケーリングとπ_intのレンジ拡張
        let validM = history_M.filter(v => !isNaN(v));
        let maxM = validM.length > 0 ? Math.max(...validM, 0.5) : 0.5;
        let scaleM = maxM * 1.2; 
        
        let scalePi_min = 0.0; // 0まで下がる可能性を考慮
        let scalePi_max = 3.0; // 少し余裕を持たせる

        // Phase boundary (λ₀ = 1.0)
        let normLambda = (1.0 - scalePi_min) / (scalePi_max - scalePi_min);
        let pxLambda = padX + normLambda * (w - padX - 10);
        this.loopCtx.beginPath();
        this.loopCtx.strokeStyle = "rgba(255, 85, 85, 0.4)";
        this.loopCtx.lineWidth = 1.5;
        this.loopCtx.setLineDash([5, 5]);
        this.loopCtx.moveTo(pxLambda, 10);
        this.loopCtx.lineTo(pxLambda, h - padY);
        this.loopCtx.stroke();
        this.loopCtx.setLineDash([]);

        this.loopCtx.font = "10px sans-serif";
        this.loopCtx.fillStyle = "rgba(255, 85, 85, 0.8)";
        this.loopCtx.fillText("λ₀", pxLambda + 5, 20);

        let lapFrames = 2500;
        let splitIdx = Math.max(0, history_pi.length - lapFrames);

        // Past loop (dotted) 過去の軌跡
        this.loopCtx.beginPath();
        this.loopCtx.strokeStyle = "rgba(74, 144, 226, 0.3)";
        this.loopCtx.lineWidth = 1.5;
        this.loopCtx.setLineDash([4, 4]);
        let firstDotted = true;
        for (let i = 0; i <= splitIdx; i++) {
            if (isNaN(history_M[i])) continue;
            let normPi = (history_pi[i] - scalePi_min) / (scalePi_max - scalePi_min);
            let px = padX + normPi * (w - padX - 10);
            let py = (h - padY) - (history_M[i] / scaleM) * (h - padY - 10);
            if (firstDotted) { this.loopCtx.moveTo(px, py); firstDotted = false; }
            else this.loopCtx.lineTo(px, py);
        }
        this.loopCtx.stroke();

        // Current loop (solid) 直近の軌跡
        this.loopCtx.beginPath();
        this.loopCtx.strokeStyle = "rgba(74, 144, 226, 1.0)";
        this.loopCtx.lineWidth = 2.5;
        this.loopCtx.setLineDash([]);
        let firstSolid = true;
        for (let i = splitIdx; i < history_pi.length; i++) {
            if (isNaN(history_M[i])) continue;
            let normPi = (history_pi[i] - scalePi_min) / (scalePi_max - scalePi_min);
            let px = padX + normPi * (w - padX - 10);
            let py = (h - padY) - (history_M[i] / scaleM) * (h - padY - 10);
            if (firstSolid) { this.loopCtx.moveTo(px, py); firstSolid = false; }
            else this.loopCtx.lineTo(px, py);
        }
        this.loopCtx.stroke();

        // Current position dot 現在地
        if (history_pi.length > 0) {
            let lastIdx = history_pi.length - 1;
            let normPi = (history_pi[lastIdx] - scalePi_min) / (scalePi_max - scalePi_min);
            let px = padX + normPi * (w - padX - 10);
            let py = (h - padY) - (history_M[lastIdx] / scaleM) * (h - padY - 10);

            this.loopCtx.beginPath();
            this.loopCtx.arc(px, py, 4, 0, Math.PI * 2);
            this.loopCtx.fillStyle = "#ffc800";
            this.loopCtx.fill();
        }
    }
}