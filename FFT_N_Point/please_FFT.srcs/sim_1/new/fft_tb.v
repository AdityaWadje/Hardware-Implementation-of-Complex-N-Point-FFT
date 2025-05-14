`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 19.03.2025 23:34:32
// Design Name: 
// Module Name: fft_tb
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



//////////////////////////////////////////////////////////////////////////////////////
module fft_tb;
     // Parameters
    parameter M = 8;           // Number of FFT points (must be power of 2)
    parameter N = 8;            // Bit width for each data point (e.g., Q1.(N-1))

    reg clk;
    reg rst;
    reg signed [(M*N)-1:0] real_in;
    reg signed [(M*N)-1:0] imag_in;
    wire signed [(M*N)-1:0] real_out;
    wire signed [(M*N)-1:0] imag_out;
    wire finish;

    // Instantiate the FFT module
    fft #(.M(M), .N(N)) uut (
        .clk(clk),
        .rst(rst),
        .real_in(real_in),
        .imag_in(imag_in),
        .real_out(real_out),
        .imag_out(imag_out),
        .finish(finish)
    );

    // Clock generation
    always #5 clk = ~clk;

    initial begin
        clk = 0;
        rst = 1;
        real_in = 0;
        imag_in = 0;
        
        #20 rst = 0;
        
        // Providing real input: {1, 2, 3, 4, 5, 6, 7, 8}
        // Imaginary input is zero
        real_in = { 8'd7, 8'd6, 8'd5, 8'd4, 8'd3, 8'd2, 8'd1, 8'd0};
//        real_in = { 4'd4, 4'd3, 4'd2, 4'd1, 4'd4, 4'd3, 4'd2, 4'd1};
//        real_in = { 16'd1, 16'd1, 16'd1, 16'd1, 16'd2, 16'd2, 16'd2, 16'd2};


//        real_in = 0;
       imag_in = { 8'd1, 8'd2, 8'd3, 8'd4, 8'd5, 8'd6, 8'd7, 8'd8};

//        imag_in = { 16'd4, 16'd3, 16'd2, 16'd1, 16'd4, 16'd3, 16'd2, 16'd1};
//        imag_in = 0;
        
    end

    // Print outputs at every positive clock edge
//    always @(posedge clk) begin
//        $display("Time = %t, Finish = %b", $time, finish);
//        $display("Real Output: %h", real_out);
//        $display("Imag Output: %h", imag_out);
//        $display("------------------------------------");
//    end

    // Stop simulation after FFT completes
  always @(posedge clk) begin
    if (finish) begin
      $stop;
    end
  end
endmodule
