// =========================================================
// AUDIO COMPRESSION USING BLOCK-BASED MDCT + QUANTIZATION
// DSP Micro-project - Scilab (v2, corrected)
//
// PIPELINE (same idea as MP3 / AAC / Vorbis / Opus):
//   1. Split audio into 50%-overlapping frames (length 2N, hop N)
//   2. Multiply each frame by a sine window
//   3. MDCT  -> N frequency coefficients per frame
//   4. QUANTIZE coefficients (bigger step for high frequencies)
//   5. ENTROPY ESTIMATE of the quantized values -> bits needed
//   6. Inverse MDCT + window + overlap-add -> reconstructed audio
//
// Compression comes from step 4 + 5: after quantization most
// coefficients become 0 or small integers, which need very few
// bits when entropy coded (Huffman/arithmetic coding).
//
// Requires Scilab 6.x (uses audioread / audiowrite).
// =========================================================

clear; clc; close(winsid());

// ---------------------------------------------------------
// STEP 1: Load audio  (CHANGE THE PATHS)
// ---------------------------------------------------------
filename      = "C:\Users\vaish\Documents\DSP_Audio_Project\my_audio.wav";
output_folder = "C:\Users\vaish\Documents\DSP_Audio_Project\Reconstructed Audio\";

if ~isdir(output_folder) then
    mkdir(output_folder);
end

try
    [xin, Fs] = audioread(filename);     // Scilab 6.x: samples x channels
    x = xin(:, 1);                       // first channel, column vector
catch
    [xl, info] = loadwave(filename);     // fallback: same as your original script
    xl = xl(1, :);
    x = xl(:);
    Fs = info(3);
end
x = x / max(abs(x));                 // normalise peak to 1
N_total = length(x);
ORIG_BITS = 16;                      // assume 16-bit PCM original

// ---------------------------------------------------------
// STEP 2: MDCT setup
// ---------------------------------------------------------
N         = 1024;                    // coefficients per frame (hop size)
frame_len = 2 * N;                   // frame length = 2N samples
hop       = N;
printf("Frame length = %.1f ms, hop = %.1f ms (Fs = %d Hz)\n", ...
       1000*frame_len/Fs, 1000*hop/Fs, Fs);

n_idx = (0:frame_len-1);             // row: time index
k_idx = (0:N-1)';                    // column: frequency index

// Sine window (Princen-Bradley: w(n)^2 + w(n+N)^2 = 1)
w = sin(%pi/frame_len * (n_idx' + 0.5));          // column, 2N x 1

// MDCT basis, N x 2N, built as an OUTER PRODUCT (works on all versions)
basis = cos((%pi/N) * ((k_idx + 0.5) * (n_idx + 0.5 + N/2)));

// Scaling sqrt(2/N) on BOTH forward and inverse makes the transform
// orthonormal, so coefficient size ~ sample size and the quantization
// error in the coefficients equals the error in the time domain.
sc = sqrt(2/N);

// Pad: N zeros in front, zeros at the end so the last frame fits
F = ceil(N_total / N) + 1;                        // number of frames
x_padded = [zeros(N,1); x; zeros((F+1)*N - N - N_total, 1)];

// ---------------------------------------------------------
// STEP 3: Forward MDCT of every frame  ->  Xall is N x F
// ---------------------------------------------------------
Xall = zeros(N, F);
for f = 1:F
    s = (f-1)*hop + 1;
    Xall(:, f) = sc * basis * (x_padded(s : s+frame_len-1) .* w);
end

// ---------------------------------------------------------
// STEP 4: Quantization settings
// Bands roughly follow how hearing works: narrow at low
// frequencies, wide at high frequencies.
// Higher bands get a LARGER step (coarser) because the ear is less
// sensitive to detail there. Set all weights to 1 for plain uniform.
// ---------------------------------------------------------
band_edges   = round(N * [0 0.04 0.10 0.20 0.35 0.60 1.0]);
band_weights = [1 1 1.5 2 3 4];
NUM_BANDS    = length(band_weights);

DEADZONE = 0.15;       // widens the zero bin -> more zeros, fewer bits

// Sweep of base step sizes (bigger step = more compression, more error)
steps    = [0.0005 0.001 0.002 0.004 0.008 0.016];
SELECTED = 4;          // which step to save and play (index into 'steps')

snr_db = zeros(1, length(steps));
bps    = zeros(1, length(steps));      // estimated bits per sample
cr     = zeros(1, length(steps));      // compression ratio
x_sel  = [];
q_sel  = [];

// ---------------------------------------------------------
// STEP 5: For every step size: quantize -> entropy -> IMDCT
// ---------------------------------------------------------
for i = 1:length(steps)

    // per-coefficient step size (column), expanded across all frames
    stepvec = zeros(N, 1);
    for b = 1:NUM_BANDS
        stepvec(band_edges(b)+1 : band_edges(b+1)) = steps(i) * band_weights(b);
    end
    stepmat = stepvec * ones(1, F);

    // Dead-zone uniform quantizer -> integers
    Q = sign(Xall) .* floor(abs(Xall) ./ stepmat + 0.5 - DEADZONE);

    // ---- Entropy estimate (ideal entropy coder, one table per band) ----
    total_bits = 0;
    for b = 1:NUM_BANDS
        qb = Q(band_edges(b)+1 : band_edges(b+1), :);
        qb = qb(:);
        t  = tabul(qb, "i");            // [value, count]
        p  = t(:,2) / sum(t(:,2));
        H  = -sum(p .* log2(p));        // bits per coefficient in this band
        total_bits = total_bits + H * length(qb);
    end
    bps(i) = total_bits / N_total;
    cr(i)  = ORIG_BITS / bps(i);

    // ---- Dequantize + inverse MDCT + overlap-add ----
    Xhat = Q .* stepmat;
    Y    = sc * (basis' * Xhat);        // 2N x F
    y_padded = zeros((F+1)*N, 1);
    for f = 1:F
        s = (f-1)*hop + 1;
        y_padded(s : s+frame_len-1) = y_padded(s : s+frame_len-1) + Y(:, f) .* w;
    end
    xr = y_padded(N+1 : N+N_total);

    err = x - xr;
    snr_db(i) = 10 * log10(sum(x.^2) / sum(err.^2));

    printf("step=%.4f | SNR=%6.2f dB | %.2f bits/sample | ratio ~ %.1f : 1 | zeros=%.1f %%\n", ...
           steps(i), snr_db(i), bps(i), cr(i), 100*sum(Q==0)/length(Q));

    if i == SELECTED then
        x_sel = xr;
        q_sel = Q;
    end
end

// ---------------------------------------------------------
// STEP 6: Measure how much was lost (same report as original script)
// ---------------------------------------------------------
x_sel = max(min(x_sel, 1), -1);          // clip only (no re-normalising)
error_signal = x - x_sel;
mse = mean(error_signal.^2);

total_coeffs = length(q_sel);
kept_coeffs  = sum(q_sel <> 0);
data_retained = (kept_coeffs / total_coeffs) * 100;

printf("Total MDCT coefficients (all frames) : %d\n", total_coeffs);
printf("Coefficients kept                    : %d\n", kept_coeffs);
printf("Data retained                        : %.2f %%\n", data_retained);
printf("Mean Squared Error                   : %e\n", mse);
printf("SNR                                  : %.2f dB\n", snr_db(SELECTED));
printf("Estimated bits per sample            : %.2f (original: %d)\n", bps(SELECTED), ORIG_BITS);
printf("Estimated compression ratio          : %.1f : 1\n", cr(SELECTED));

// ---------------------------------------------------------
// STEP 7: Save and listen
// ---------------------------------------------------------
output_filename = output_folder + "reconstructed_audio_mdct.wav";
try
    audiowrite(output_filename, x_sel, Fs);
catch
    savewave(output_filename, x_sel', Fs);
end
printf("\nSaved reconstructed audio to: %s\n", output_filename);

printf("Playing RECONSTRUCTED (MDCT-compressed) audio...\n");
playsnd(x_sel', Fs);

// ---------------------------------------------------------
// STEP 8: Plots (same 3 plots as original script)
// ---------------------------------------------------------
scf(1);
subplot(3,1,1);
plot(x, 'b');
title('Original Audio Signal');
xlabel('Sample'); ylabel('Amplitude');

subplot(3,1,2);
plot(x_sel, 'r');
title('Reconstructed Audio (Block MDCT, quantized, step = ' + string(steps(SELECTED)) + ')');
xlabel('Sample'); ylabel('Amplitude');

subplot(3,1,3);
plot(error_signal, 'k');
title('Error (Original - Reconstructed)');
xlabel('Sample'); ylabel('Amplitude');

// Extra: rate-distortion curve from the step-size sweep
scf(2);
plot(bps, snr_db, 'o-');
title('Rate-Distortion: quality vs. bits per sample');
xlabel('Estimated bits per sample'); ylabel('SNR (dB)'); xgrid();
