`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 

// Design Name: 
// Module Name: fft4
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: 
// 
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////
module dft_4point (
    input [31:0] real_in,  // 32-bit input for real parts (4 numbers, each 8 bits)
    input [31:0] imag_in,  // 32-bit input for imaginary parts (4 numbers, each 8 bits)
    output reg [31:0] real_out, // 32-bit output for real parts (4 numbers, each 8 bits)
    output reg [31:0] imag_out  // 32-bit output for imaginary parts (4 numbers, each 8 bits)
);

    // Twiddle factor coefficients (real and imaginary parts as fixed-point scaled by 2^14)
    reg signed [15:0] W_REAL [0:3]; // cos(0), cos(90), cos(180), cos(270)
    reg signed [15:0] W_IMAG [0:3]; // sin(0), sin(90), sin(180), sin(270)

    // Input signals unpacked into arrays
    reg signed [7:0] real_in_arr [3:0];
    reg signed [7:0] imag_in_arr [3:0];

    // Output signals packed into arrays
    reg signed [15:0] real_out_arr [3:0];
    reg signed [15:0] imag_out_arr [3:0];

    // Temporary accumulators for outputs
    reg signed [31:0] temp_real, temp_imag;
    integer k, n, twiddle_index; // Loop indices

    initial begin
        // Initialize twiddle factors
        W_REAL[0] = 16384;  W_IMAG[0] = 0;       // cos(0), sin(0)
        W_REAL[1] = 0;      W_IMAG[1] = -16384;  // cos(90), sin(90)
        W_REAL[2] = -16384; W_IMAG[2] = 0;       // cos(180), sin(180)
        W_REAL[3] = 0;      W_IMAG[3] = 16384;   // cos(270), sin(270)
    end

    always @(*) begin
        // Unpack 32-bit inputs into 4x8-bit arrays
        real_in_arr[0] = real_in[31:24];
        real_in_arr[1] = real_in[23:16];
        real_in_arr[2] = real_in[15:8];
        real_in_arr[3] = real_in[7:0];

        imag_in_arr[0] = imag_in[31:24];
        imag_in_arr[1] = imag_in[23:16];
        imag_in_arr[2] = imag_in[15:8];
        imag_in_arr[3] = imag_in[7:0];

        // Initialize output arrays
        for (k = 0; k < 4; k = k + 1) begin
            temp_real = 0;
            temp_imag = 0;

            // Perform DFT calculations for each output (k)
            for (n = 0; n < 4; n = n + 1) begin
            twiddle_index = (k * n) % 4;

                // Real part computation
                temp_real = temp_real +
                            (real_in_arr[n] * W_REAL[twiddle_index] - imag_in_arr[n] * W_IMAG[twiddle_index]);

                // Imaginary part computation
                temp_imag = temp_imag +
                            (real_in_arr[n] * W_IMAG[twiddle_index] + imag_in_arr[n] * W_REAL[twiddle_index]);
            end

            // Scale down by 2^14 and assign to output arrays
            real_out_arr[k] = temp_real >>> 14;
            imag_out_arr[k] = temp_imag >>> 14;
        end

        // Pack output arrays into 32-bit outputs
        real_out = {real_out_arr[0][7:0], real_out_arr[1][7:0], real_out_arr[2][7:0], real_out_arr[3][7:0]};
        imag_out = {imag_out_arr[0][7:0], imag_out_arr[1][7:0], imag_out_arr[2][7:0], imag_out_arr[3][7:0]};
    end

endmodule
