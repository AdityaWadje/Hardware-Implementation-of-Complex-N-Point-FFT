`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 19.03.2025 21:19:11
// Design Name: 
// Module Name: fft
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


module fft #(
    parameter M = 8,           // Number of FFT points (must be power of 2)
    parameter N = 8           // Bit width for each data point (e.g., Q1.(N-1))
)(
    input clk,
    input rst,
//    input signed [(M*N)-1:0] real_in,
//    input signed [(M*N)-1:0] imag_in,
//    output reg signed [(M*N)-1:0] real_out,
//    output reg signed [(M*N)-1:0] imag_out,
    output reg signed [(N)-1:0] real_out,
    output reg signed [(N)-1:0] imag_out,
    output reg finish // Finish flag
);
    // Unpacking input signals into arrays
    reg signed [N-1:0] real_data [0:M-1];
    reg signed [N-1:0] imag_data [0:M-1];
    reg signed [N-1:0] real_data_out [0:M-1];
    reg signed [N-1:0] imag_data_out [0:M-1];
    reg finish1;   
    wire signed [N-1:0] real_part_u [0:(M/2*$clog2(M))-1];
    wire signed [N-1:0] img_part_u[0:(M/2*$clog2(M))-1];
    wire signed [N-1:0] real_part_l [0:(M/2*$clog2(M))-1];
    wire signed [N-1:0] img_part_l[0:(M/2*$clog2(M))-1];
    reg signed [N-1:0] buffer_real [0:(M*$clog2(M))-1];
    reg signed [N-1:0] buffer_imag [0:(M*$clog2(M))-1];
//    reg signed [N-1:0] buffer_real_l [0:(M/2*$clog2(M))-1];
//    reg signed [N-1:0] buffer_imag_l [0:(M/2*$clog2(M))-1];
    // Processing FFT using Optimized Radix-2 Butterfly Computation
    localparam log_M = $clog2(M);
    localparam log_R = $clog2(log_M*M/2);
    reg [($clog2(M)*M/2):0] finish_radix ;
    reg finish_seperation;
    reg [log_R:0] finish_count;
    reg [2:0] state;
    reg [M/2*$clog2(M/2): 0]rst_radix;
    reg [M:0]count;
    localparam  RST_ON = 3'b000, RST_OFF = 3'b001, EXCHANGE = 3'b010, READY = 3'b011, DONE = 3'b100, OUT_READY = 3'b101, DON2 = 3'b110,  COMBINE = 3'b111 ;
//    reg real_in = {8'd4, 8'd3, 8'd2, 8'd1, 8'd4, 8'd3, 8'd2, 8'd1};
//    reg imag_in = {8'd4, 8'd3, 8'd2, 8'd1, 8'd4, 8'd3, 8'd2, 8'd1};


    integer j, t, y, rev, temp;
    always @(posedge clk or posedge rst) begin
        if (rst) begin
            for (j = 0; j < M; j = j + 1) begin
                real_data[j] <= 0;
                imag_data[j] <= 0;
                real_data_out[j] <= 0;
                imag_data_out[j] <= 0;              
            end
            for (t = 1; t < ($clog2(M)*M/2); t = t + 1) begin
                finish_radix[t] <= 1'b0;
                rst_radix[t] <= 1'b0;
//                $display("finish_radix[%0d]: %d", t, finish_radix[t]);
            end  
            finish <= 0; // Reset finish flag
            finish_seperation <= 0;     
            finish1 <= 0;
            finish_radix[0]  = 1'b1;  // Assign binary 1 to the register
//            $display("finish_radix[%0d]: %d", 0, finish_radix[0]);
            finish_count <= 1;
            rst_radix[0] <= 1; 
            state <= RST_ON;
            real_out <= 0;
            imag_out <= 0;
            count <= 0; 

        end else if (!finish_seperation) begin
            for (j = 0; j < M; j = j + 1) begin
//                real_data[j] <= real_in[(j*N) +: N];
//                imag_data[j] <= imag_in[(j*N) +: N];      
                real_data[j] <= j;
                imag_data[j] <= M-j; 
            end
            finish_seperation <= 1;
            finish <= 0; // Reset at start
        end
    end


reg [$clog2(M)-1:0]out[0:M-1];
//reg finish_radix_bit;
wire radix_finish;

  
  genvar s, g, i, r, r_g;
//  genvar butterfly_index;
  genvar twiddle_multiplier;
  genvar j_bits, g_bits;
  
  // Maximum size for the temporary arrays:
  // For inner indices: max size is 2^(M-1) (when s = 0)
  // For group indices: max size is N (when L=1, though typically groups are fewer)
  reg ordered_i [0:(2**(log_M-1))-1];
  reg ordered_g [0:log_M-1];

  // Bit reversal function: reverses 'bits' bits of x.
  function integer bit_reverse;
    input integer x;
    input integer bits;
    begin
      bit_reverse = 0;
      for (j = 0; j < bits; j = j + 1)
        if (x & (2**j))
          bit_reverse = bit_reverse | (1 << (bits - 1 - j));
    end
  endfunction

generate
    for (s = 0; s < log_M; s = s + 1) begin : stage_loop
      localparam L = (2**(log_M - s));           // Group length
      localparam num_groups = M / L;           // Number of groups
      localparam num_butterflies = L / 2;      // Butterflies per group
      localparam twiddle_multiplier = 2**(log_M - s - 1);

      // Generate butterflies inside each group
      for (i = 0; i < num_butterflies; i = i + 1) begin : butterfly_loop
        for (g = 0; g < num_groups; g = g + 1) begin : group_loop
          wire [log_M-1:0] top, bottom;
          wire [4:0] twiddle;

          assign top = (g * L) + i;
          assign bottom = top + num_butterflies;
          assign twiddle = g * twiddle_multiplier;
          
//          assign twiddle = g * twiddle_multiplier*(16/M);
//          initial begin
//            $display("%0d-%0d-%0d-%0d", top, bottom, twiddle, (s * num_butterflies * num_groups) + (j * num_groups) + g);
            radix2_butterfly #(.M(M), .N(N)) uut (
                    .a_real(real_data[top]),
                    .a_imag(imag_data[top]), 
                    .b_real(real_data[bottom]), 
                    .b_imag(imag_data[bottom]),
                    .out1_real(real_part_u[(s * num_butterflies * num_groups) + (i * num_groups) + g]), 
                    .out1_imag(img_part_u[(s * num_butterflies * num_groups) + (i * num_groups) + g]),
                    .out2_real(real_part_l[(s * num_butterflies * num_groups) + (i * num_groups) + g]), 
                    .out2_imag(img_part_l[(s * num_butterflies * num_groups) + (i * num_groups) + g]),.fft_index(twiddle),
                    .done(radix_finish), .clk(clk), .reset(rst_radix[finish_count-1]),.start(1)

                );
//          end
        end
      end
    end

  endgenerate

            always @(posedge clk) begin
                finish_radix[finish_count] = radix_finish;
                case (state)
                    RST_ON: begin
                        rst_radix[finish_count] <= 1;
                        state <= RST_OFF;
//                        $display("state: %d", state);
                    end 
                    RST_OFF: begin
                        rst_radix[finish_count - 1] <= 0;
                        if (finish_seperation && finish_radix[finish_count]) begin
                            state <= EXCHANGE;
                        end
                    end    
                    
                    EXCHANGE: begin
                    for (t = 0; t < log_M; t = t + 1) begin
                        for (j = 0; j < (2**(log_M - t)) / 2; j = j + 1) begin
                            for (y = 0; y < M /(2**(log_M - t)); y = y + 1) begin
                                buffer_real[(y * (2**(log_M - t))) + j + M*t] = real_part_u[(t * ((2**(log_M - t)) / 2) * (M /(2**(log_M - t)))) + (j * (M /(2**(log_M - t)))) + y];
                                buffer_imag[(y * (2**(log_M - t))) + j + M*t] = img_part_u[(t * ((2**(log_M - t)) / 2) * (M /(2**(log_M - t)))) + (j * (M /(2**(log_M - t)))) + y];
                                buffer_real[(y * (2**(log_M - t))) + j + M*t + (2**(log_M - t)) / 2] = real_part_l[(t * ((2**(log_M - t)) / 2) * (M /(2**(log_M - t)))) + (j * (M /(2**(log_M - t)))) + y];
                                buffer_imag[(y * (2**(log_M - t))) + j + M*t + (2**(log_M - t)) / 2] = img_part_l[(t * ((2**(log_M - t)) / 2) * (M /(2**(log_M - t)))) + (j * (M /(2**(log_M - t)))) + y];
                                $display("t = %0d, j = %0d, y = %0d", t, j, y);

                                $display("buffer_real[%0d] = %d, buffer_img[%0d] = %d", (y * (2**(log_M - t))) + j + M*t, buffer_real[(y * (2**(log_M - t))) + j + M*t],(y * (2**(log_M - t))) + j + M*t, buffer_imag[(y * (2**(log_M - t))) + j + M*t]);
                                $display("buffer_real[%0d] = %d, buffer_img[%0d] = %d", (y * (2**(log_M - t))) + j + M*t + (2**(log_M - t)) / 2, buffer_real[(y * (2**(log_M - t))) + j + M*t + (2**(log_M - t)) / 2], (y * (2**(log_M - t))) + j + M*t + (2**(log_M - t)) / 2, buffer_imag[(y * (2**(log_M - t))) + j + M*t + (2**(log_M - t)) / 2]); 
                            end
                        end
                    end
                        state <= READY;
                        $display("Time2: %0t ns, state: %0d", $time, state);

//                        $display("Real[%0d] = %d, Imag[%0d] = %d", p+i, real_data[p+i], p+i, imag_data[p+i]);
//                        $display("Real[%0d] = %d, Imag[%0d] = %d", p+i+(step / 2), real_data[p+i+(step / 2)], p+i, imag_data[p+i+(step / 2)]);
                    end
                    READY: begin
                        finish_count = finish_count + 1;
                        for (t = 0; t < M; t = t + 1) begin
                                real_data [t] <= buffer_real[(finish_count - 2)*M + t];
                                imag_data [t] <= buffer_imag[(finish_count - 2)*M + t];

                        end
                       if (finish_count == $clog2(M)+1)  begin
                       

                        state <= DONE;   
                        finish1 <= 1;
                        end  else begin
                                state <= RST_ON;
                        end 
                    end
                    DONE: begin
                        state <= OUT_READY;
                        for (t = 0; t < M; t = t + 1) begin                        
                            $display("real_data[%0d] = %d, img_data[%0d] = %d", t, real_data[t], t, imag_data[t]);
                        end
                        $display("Time3: %0t ns, finish_count: %d, state: %0d", $time, finish_count, state);  

                    end
                    OUT_READY: begin
                    for (t = 0; t < M; t = t + 1) begin : gen_block
                    rev = 0;
                    temp = t;
                    for (j=0; j< $clog2(M); j = j+1) begin
                        rev = (rev << 1) | (temp & 1);  // Shift left and take LSB
                        temp = temp >> 1;  // Shift right
                    end
                    if ( finish1) begin
                            real_data_out[t] <= real_data[rev];
                            imag_data_out[t]  <= imag_data[rev];
                            count <= count+1;

                    end 
                    end                          
                    state <= DON2;

                    end

                    DON2: begin
                        if (count == M) begin
                            finish <= 1;  // Set finish signal when all iterations are complete
                            finish1 <= 0;  // Set finish signal when all iterations are complete
                            state <= COMBINE;
                            real_out <= real_data_out[count-1];
                            imag_out <= imag_data_out[count-1];
                        end else begin
                            state <= OUT_READY;
                            real_out <= real_data_out[count-1];
                            imag_out <= imag_data_out[count-1];
                        end
                    $display("Time3: %0t ns, count: %d, state: %0d, real: %0d, imag: %0d", $time, count, state, real_data_out[count-1], imag_data_out[count-1]);  
                    end 
                    COMBINE: begin
                        if (finish) begin
                              $display("Time: %0t ns, finish: %d, state: %0d", $time, finish, state);  
                        end
                    end
               endcase
            end

endmodule

module radix2_butterfly #(
    parameter M = 8,           // Number of FFT points (must be power of 2)
    parameter N = 8           // Bit width for each data point (e.g., Q1.(N-1))
)(
    input clk,
    input reset,
    input start,
    input [$clog2(M)-1:0] fft_index,  // FFT index for the twiddle factor
    // Complex inputs
    input  signed [N-1:0] a_real,
    input  signed [N-1:0] a_imag,
    input  signed [N-1:0] b_real,
    input  signed [N-1:0] b_imag,
    // Butterfly outputs
    output reg signed [N-1:0] out1_real,
    output reg signed [N-1:0] out1_imag,
    output reg signed [N-1:0] out2_real,
    output reg signed [N-1:0] out2_imag,
    output reg done
);

    // Signals from the CORDIC FFT twiddle generator
    wire signed [N-1:0] tw_cos;
    wire signed [N-1:0] tw_sin;
    wire cordic_done;

    // Instantiate the updated CORDIC FFT twiddle generator.
    cordic_dynamic #(
        .M(M),
        .N(N)
    ) cordic_inst (
        .clk(clk),
        .reset(reset),
        .start(start),
        .addr(fft_index),
        .cos_out(tw_cos),
        .sin_out(tw_sin),
        .done(cordic_done)
    );

    // Intermediate registers for the twiddle multiplication result.
    // Multiplication results in double width.
    reg signed [2*N-1:0] prod_real, prod_imag;
    reg signed [2*N-1:0] b_twiddled_real, b_twiddled_imag;
    localparam SCALE_FACTOR = N - 1; // For Q1.15 fixed-point format

    // A simple state machine for controlling the butterfly computation.
    reg [1:0] state;
    localparam IDLE       = 2'b00,
               WAIT_CORDIC = 2'b01,
               COMPUTE    = 2'b10,
               OUTPUT     = 2'b11;

    always @(posedge clk or posedge reset) begin
        if (reset) begin
            state       <= IDLE;
            done        <= 0;
            out1_real   <= 0;
            out1_imag   <= 0;
            out2_real   <= 0;
            out2_imag   <= 0;
        end else begin
            case(state)
                IDLE: begin
                    done <= 0;
                    if (start) begin
                        state <= WAIT_CORDIC;
                    end
                end

                // Wait until the CORDIC module has generated the twiddle factor.
                WAIT_CORDIC: begin
                    if (cordic_done) begin
                        state <= COMPUTE;
                    end
                end

                // Compute the complex multiplication: b * twiddle.
                COMPUTE: begin
                    // Complex multiplication:
                    //   Real part = b_real*tw_cos - b_imag*tw_sin
                    //   Imag part = b_real*tw_sin + b_imag*tw_cos
                    prod_real = b_real * tw_cos - b_imag * tw_sin;
                    prod_imag = b_real * tw_sin + b_imag * tw_cos;

                    // Adjust for fixed-point scaling using rounding.
                    b_twiddled_real = (prod_real >= 0) ? 
                        ((prod_real + (1 << (SCALE_FACTOR - 1))) >>> SCALE_FACTOR) :
                        -(( -prod_real + (1 << (SCALE_FACTOR - 1))) >>> SCALE_FACTOR);

                    b_twiddled_imag = (prod_imag >= 0) ? 
                        ((prod_imag + (1 << (SCALE_FACTOR - 1))) >>> SCALE_FACTOR) :
                        -(( -prod_imag + (1 << (SCALE_FACTOR - 1))) >>> SCALE_FACTOR);
                    
                    out1_real <= a_real + b_twiddled_real;
                    out1_imag <= a_imag + b_twiddled_imag;
                    out2_real <= a_real - b_twiddled_real;
                    out2_imag <= a_imag - b_twiddled_imag;

                    state <= OUTPUT;
                end

                // Compute the butterfly addition/subtraction.
                OUTPUT: begin
                    // Butterfly outputs:
                    // out1 = A + (B * twiddle)
                    // out2 = A - (B * twiddle)
//                    out1_real <= a_real + b_twiddled_real;
//                    out1_imag <= a_imag + b_twiddled_imag;
//                    out2_real <= a_real - b_twiddled_real;
//                    out2_imag <= a_imag - b_twiddled_imag;
                    
                    $display("Time7: %0t ns", $time);
                    $display("rReal1 = %d, rImag1 = %d", a_real, a_imag);
                    $display("rReal2 = %d, rImag2 = %d", b_real, b_imag);
                    $display("rReal_1 = %d, rImag_1 = %d", out1_real, out1_imag);
                    $display("rReal_2 = %d, rImag_2 = %d", out2_real, out2_imag);
                    $display("wReal = %d, wImag = %d, index = %0d", tw_cos, tw_sin, fft_index);  
                    done <= 1;
                    state <= IDLE;
                end

                default: state <= IDLE;
            endcase
        end
    end

endmodule



module cordic_dynamic #(
    parameter N = 8,           
    parameter M =8,
    parameter ADDR_WIDTH = $clog2(M)  // 16-bit angle representation
)(
    input clk,
    input reset,
    input start,
    input  wire [ADDR_WIDTH-1:0] addr,
    output reg signed [N-1:0] cos_out,
    output reg signed [N-1:0] sin_out,
    output reg done
);

    // Scaling factor based on fixed-point representation, e.g., Q1.(DATA_WIDTH-1)
  localparam integer SCALE = (1 << (N-1)) - 1;

  // ROM arrays for cosine and sine values.
  reg signed [N-1:0] cos_rom [0:M-1];
  reg signed [N-1:0] sin_rom [0:M-1];

  // Loop variable and temporary real value for computation.
  integer i;
  real angle;

  // Precompute the twiddle factors:
  // For each index i, calculate the angle -2*pi*i/N and convert the cosine and sine values
  // into fixed-point numbers using the scaling factor.
  initial begin
    for (i = 0; i < M; i = i + 1) begin
      angle = -2.0 * 3.14159265358979323846 * i / M;
      cos_rom[i] = $rtoi($cos(angle) * SCALE);
      sin_rom[i] = $rtoi($sin(angle) * SCALE);
    end
  end

  // On each clock cycle, output the corresponding twiddle factor based on the address.
  always @(posedge clk) begin
    cos_out <= cos_rom[addr];
    sin_out <= sin_rom[addr];
    done   <= 1;

  end
endmodule
