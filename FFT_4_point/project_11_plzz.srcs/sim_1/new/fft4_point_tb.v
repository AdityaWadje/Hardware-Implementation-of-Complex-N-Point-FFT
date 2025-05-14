`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 

// Design Name: 
// Module Name: fft4_point_tb
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
module dft_4point_tb;

    reg [31:0] real_in;  // 32-bit input for real parts
    reg [31:0] imag_in;  // 32-bit input for imaginary parts
    wire [31:0] real_out; // 32-bit output for real parts
    wire [31:0] imag_out; // 32-bit output for imaginary parts

    // Instantiate the DFT module
    dft_4point uut (
        .real_in(real_in),
        .imag_in(imag_in),
        .real_out(real_out),
        .imag_out(imag_out)
    );

    initial begin
        // Test case 1: Example inputs
        real_in = {8'd9, 8'd10, 8'd11, 8'd12}; // Real parts: 9, 10, 11, 12
        imag_in = {8'd0, 8'd0, 8'd0, 8'd0};   // Imaginary parts: 0, 0, 0, 0
        #1; // Wait for computation
        $display("Test 1: Real Out: %h, Imag Out: %h", real_out, imag_out);

        // Test case 2: Complex inputs
        real_in = {8'd1, 8'd2, 8'd3, 8'd4};  // Real parts: 1, 2, 3, 4
        imag_in = {8'd4, 8'd3, 8'd2, 8'd1};  // Imaginary parts: 4, 3, 2, 1
        #1; // Wait for computation
        $display("Test 2: Real Out: %h, Imag Out: %h", real_out, imag_out);

        // Test case 3: All zeros
        real_in = {8'd0, 8'd0, 8'd0, 8'd0};  // Real parts: 0, 0, 0, 0
        imag_in = {8'd0, 8'd0, 8'd0, 8'd0};  // Imaginary parts: 0, 0, 0, 0
        #1; // Wait for computation
        $display("Test 3: Real Out: %h, Imag Out: %h", real_out, imag_out);

        $finish;
    end

endmodule
