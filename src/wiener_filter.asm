# =============================================================================
# Wiener filter (filtering & MMSE) in MIPS32 assembly - runs on MARS 4.5
#
# Reads   input.txt   (noisy signal x[n], 10 values, 1 decimal)
#         desired.txt (desired signal d[n], 10 values, 1 decimal)
# Writes  output.txt  and prints the same two lines to the console:
#         Filtered output: y[0] ... y[9]
#         MMSE: <value>
#
# Pipeline: parse ASCII -> double | correlation (GXX, GDX) | Toeplitz R
#           | Gaussian elimination + back substitution (R h = GDX)
#           | convolution y = h * x | MMSE
#
# File paths are relative to MARS's working directory (the folder MARS is
# launched from). Use absolute paths if your setup needs them.
# =============================================================================

.data
    .align 3
    N_double:    .double 10.0
    zero_double: .double 0.0
    ten_double:  .double 10.0

    filename_input:   .asciiz "input.txt"
    filename_desired: .asciiz "desired.txt"
    filename_output:  .asciiz "output.txt"

    msgFileErr: .asciiz "Error : Cannot Open File !"
    msgSizeErr: .asciiz "Error : size not match "

    msgY:       .asciiz "Filtered output: "
    msgMMSE:    .asciiz "\nMMSE: "
    space:      .asciiz " "
    dot:        .asciiz "."
    minus:      .asciiz "-"

    .align 3
    buffer: .space 1024

    .align 3
    input_signal:         .space 80
    desired_signal:       .space 80
    output_signal:        .space 80
    R:                    .space 800
    optimize_coefficient: .space 80
    GXX:                  .space 80
    GDX:                  .space 80

    .align 2
    size_input_signal:   .word 0
    size_desired_signal: .word 0

    one_char:   .space 1
    int_buffer: .space 20

.text
.globl main

main:
    la $a0, filename_input
    la $a1, input_signal
    la $a2, size_input_signal
    jal loadInputFile
    beqz $v0, file_error

    la $a0, filename_desired
    la $a1, desired_signal
    la $a2, size_desired_signal
    jal loadInputFile
    beqz $v0, file_error

    lw $t0, size_input_signal
    lw $t1, size_desired_signal

    li $t2, 10
    bne $t0, $t2, size_error
    bne $t1, $t2, size_error
    bne $t0, $t1, size_error

    li $s7, 0

main_corr_loop:
    bge $s7, 10, main_after_corr
    move $a0, $s7
    jal build_Correlation
    addi $s7, $s7, 1
    j main_corr_loop

main_after_corr:
    jal build_R_matrix
    jal calculate_H_matrix
    jal calculate_output
    jal calculate_MMSE

    li $v0, 13
    la $a0, filename_output
    li $a1, 1
    li $a2, 0
    syscall
    bltz $v0, file_error
    move $s6, $v0

    la $a0, msgY
    jal write_string_file

    li $v0, 4
    la $a0, msgY
    syscall

    li $s5, 0

write_Y_loop:
    bge $s5, 10, write_MMSE_start

    sll $t1, $s5, 3
    la $t2, output_signal
    add $t2, $t2, $t1
    l.d $f12, 0($t2)

    mov.d $f30, $f12
    jal print_double_1digit_console

    mov.d $f12, $f30
    jal write_double_1digit_file

    la $a0, space
    jal write_string_file

    li $v0, 4
    la $a0, space
    syscall

    addi $s5, $s5, 1
    j write_Y_loop

write_MMSE_start:
    la $a0, msgMMSE
    jal write_string_file

    li $v0, 4
    la $a0, msgMMSE
    syscall

    mov.d $f12, $f0
    jal print_double_1digit_console

    mov.d $f12, $f0
    jal write_double_1digit_file

    li $v0, 16
    move $a0, $s6
    syscall

    li $v0, 10
    syscall

file_error:
    li $v0, 4
    la $a0, msgFileErr
    syscall
    li $v0, 10
    syscall

size_error:
    li $v0, 4
    la $a0, msgSizeErr
    syscall
    li $v0, 10
    syscall

# ----------------------------------------------------------------------------------
# -----------------------------< Load Input File >----------------------------------
# ----------------------------------------------------------------------------------

loadInputFile:
    addi $sp, $sp, -4
    sw $ra, 0($sp)

    move $s1, $a1
    move $s2, $a2

    li $v0, 13
    li $a1, 0
    li $a2, 0
    syscall

    bltz $v0, loadInputFile_error
    move $s0, $v0

    li $t0, 0
    sw $t0, 0($s2)

    li $v0, 14
    move $a0, $s0
    la $a1, buffer
    li $a2, 1024
    syscall

    move $s3, $v0

    li $v0, 16
    move $a0, $s0
    syscall

    li $t0, 0
    li $t1, 0
    li $t2, 0
    li $t3, 1
    li $t4, 0
    li $t5, 0
    li $s4, 1

parse_loop:
    bge $t0, $s3, parse_end

    la $t6, buffer
    add $t6, $t6, $t0
    lb $t7, 0($t6)

    beq $t7, 32, save_number_if_needed
    beq $t7, 10, save_number_if_needed
    beq $t7, 13, save_number_if_needed
    beq $t7, 46, mark_after_dot
    beq $t7, 45, mark_negative

    blt $t7, 48, parse_next_char
    bgt $t7, 57, parse_next_char

    addi $t7, $t7, -48
    li $t4, 1

    beqz $t5, add_int_part

add_frac_part:
    mul $t2, $t2, 10
    add $t2, $t2, $t7
    mul $t3, $t3, 10
    j parse_next_char

add_int_part:
    mul $t1, $t1, 10
    add $t1, $t1, $t7
    j parse_next_char

mark_after_dot:
    li $t5, 1
    li $t4, 1
    j parse_next_char

mark_negative:
    li $s4, -1
    li $t4, 1
    j parse_next_char

save_number_if_needed:
    beqz $t4, parse_next_char
    jal save_parsed_number

    li $t1, 0
    li $t2, 0
    li $t3, 1
    li $t4, 0
    li $t5, 0
    li $s4, 1

parse_next_char:
    addi $t0, $t0, 1
    j parse_loop

parse_end:
    beqz $t4, loadInputFile_success
    jal save_parsed_number

loadInputFile_success:
    li $v0, 1
    lw $ra, 0($sp)
    addi $sp, $sp, 4
    jr $ra

loadInputFile_error:
    li $v0, 0
    lw $ra, 0($sp)
    addi $sp, $sp, 4
    jr $ra

save_parsed_number:
    mtc1 $t1, $f0
    cvt.d.w $f0, $f0

    mtc1 $t2, $f2
    cvt.d.w $f2, $f2

    mtc1 $t3, $f4
    cvt.d.w $f4, $f4

    div.d $f2, $f2, $f4
    add.d $f6, $f0, $f2

    mtc1 $s4, $f8
    cvt.d.w $f8, $f8
    mul.d $f6, $f6, $f8

    lw $t8, 0($s2)
    sll $t9, $t8, 3
    add $t9, $s1, $t9
    s.d $f6, 0($t9)

    addi $t8, $t8, 1
    sw $t8, 0($s2)

    jr $ra

# ----------------------------------------------------------------------------------
# ----------------------------< Build Correlation >---------------------------------
# ----------------------------------------------------------------------------------

build_Correlation:
    li $t0, 0
    l.d $f0, zero_double
    l.d $f2, zero_double

corr_loop:
    bge $t0, 10, corr_done

    sub $t1, $t0, $a0
    bltz $t1, corr_next

    sll $t2, $t0, 3
    la $t3, input_signal
    add $t3, $t3, $t2
    l.d $f4, 0($t3)

    sll $t2, $t1, 3
    la $t3, input_signal
    add $t3, $t3, $t2
    l.d $f6, 0($t3)

    sll $t2, $t0, 3
    la $t3, desired_signal
    add $t3, $t3, $t2
    l.d $f8, 0($t3)

    mul.d $f10, $f4, $f6
    add.d $f0, $f0, $f10

    mul.d $f10, $f8, $f6
    add.d $f2, $f2, $f10

corr_next:
    addi $t0, $t0, 1
    j corr_loop

corr_done:
    l.d $f14, N_double
    div.d $f0, $f0, $f14
    div.d $f2, $f2, $f14

    sll $t1, $a0, 3

    la $t2, GXX
    add $t2, $t2, $t1
    s.d $f0, 0($t2)

    la $t2, GDX
    add $t2, $t2, $t1
    s.d $f2, 0($t2)

    jr $ra

# ----------------------------------------------------------------------------------
# -----------------------------< Build R Matrix>------------------------------------
# ----------------------------------------------------------------------------------

build_R_matrix:
    li $t0, 0

R_i_loop:
    bge $t0, 10, R_done
    li $t1, 0

R_j_loop:
    bge $t1, 10, R_next_i

    sub $t2, $t0, $t1
    bgez $t2, R_abs_done
    sub $t2, $zero, $t2

R_abs_done:
    sll $t3, $t2, 3
    la $t4, GXX
    add $t4, $t4, $t3
    l.d $f0, 0($t4)

    mul $t5, $t0, 10
    add $t5, $t5, $t1
    sll $t5, $t5, 3
    la $t6, R
    add $t6, $t6, $t5
    s.d $f0, 0($t6)

    addi $t1, $t1, 1
    j R_j_loop

R_next_i:
    addi $t0, $t0, 1
    j R_i_loop

R_done:
    jr $ra

# ----------------------------------------------------------------------------------
# -----------------------< Build Optimize Coefficient>------------------------------
# ----------------------------------------------------------------------------------

calculate_H_matrix:
    li $t0, 0

gauss_k:
    bge $t0, 10, back_substitution
    addi $t1, $t0, 1

gauss_i:
    bge $t1, 10, gauss_next_k

    mul $t2, $t1, 10
    add $t2, $t2, $t0
    sll $t2, $t2, 3
    la $t3, R
    add $t3, $t3, $t2
    l.d $f0, 0($t3)

    mul $t4, $t0, 10
    add $t4, $t4, $t0
    sll $t4, $t4, 3
    la $t5, R
    add $t5, $t5, $t4
    l.d $f2, 0($t5)

    div.d $f4, $f0, $f2
    move $t6, $t0

gauss_j:
    bge $t6, 10, update_GDX

    mul $t7, $t1, 10
    add $t7, $t7, $t6
    sll $t7, $t7, 3
    la $t8, R
    add $t8, $t8, $t7
    l.d $f6, 0($t8)

    mul $t9, $t0, 10
    add $t9, $t9, $t6
    sll $t9, $t9, 3
    la $s0, R
    add $s0, $s0, $t9
    l.d $f8, 0($s0)

    mul.d $f10, $f4, $f8
    sub.d $f6, $f6, $f10
    s.d $f6, 0($t8)

    addi $t6, $t6, 1
    j gauss_j

update_GDX:
    sll $t2, $t1, 3
    la $t3, GDX
    add $t3, $t3, $t2
    l.d $f6, 0($t3)

    sll $t4, $t0, 3
    la $t5, GDX
    add $t5, $t5, $t4
    l.d $f8, 0($t5)

    mul.d $f10, $f4, $f8
    sub.d $f6, $f6, $f10
    s.d $f6, 0($t3)

    addi $t1, $t1, 1
    j gauss_i

gauss_next_k:
    addi $t0, $t0, 1
    j gauss_k

back_substitution:
    li $t0, 9

back_i:
    bltz $t0, calculate_H_done

    sll $t1, $t0, 3
    la $t2, GDX
    add $t2, $t2, $t1
    l.d $f0, 0($t2)

    addi $t3, $t0, 1

back_j:
    bge $t3, 10, save_H

    mul $t4, $t0, 10
    add $t4, $t4, $t3
    sll $t4, $t4, 3
    la $t5, R
    add $t5, $t5, $t4
    l.d $f2, 0($t5)

    sll $t6, $t3, 3
    la $t7, optimize_coefficient
    add $t7, $t7, $t6
    l.d $f4, 0($t7)

    mul.d $f6, $f2, $f4
    sub.d $f0, $f0, $f6

    addi $t3, $t3, 1
    j back_j

save_H:
    mul $t4, $t0, 10
    add $t4, $t4, $t0
    sll $t4, $t4, 3
    la $t5, R
    add $t5, $t5, $t4
    l.d $f2, 0($t5)

    div.d $f0, $f0, $f2

    sll $t6, $t0, 3
    la $t7, optimize_coefficient
    add $t7, $t7, $t6
    s.d $f0, 0($t7)

    addi $t0, $t0, -1
    j back_i

calculate_H_done:
    jr $ra

# ----------------------------------------------------------------------------------
# -----------------------< Calculate Output Signal>---------------------------------
# ----------------------------------------------------------------------------------

round_double_1digit:
    l.d $f20, ten_double
    mul.d $f12, $f12, $f20

    round.w.d $f22, $f12
    cvt.d.w $f12, $f22

    div.d $f12, $f12, $f20

    jr $ra

calculate_output:
    addi $sp, $sp, -4
    sw $ra, 0($sp)

    li $t0, 0

output_y_loop:
    bge $t0, 10, output_done

    l.d $f0, zero_double
    li $t1, 0

output_h_loop:
    bge $t1, 10, save_output_y

    sub $t2, $t0, $t1
    bltz $t2, save_output_y

    sll $t3, $t1, 3
    la $t4, optimize_coefficient
    add $t4, $t4, $t3
    l.d $f2, 0($t4)

    sll $t5, $t2, 3
    la $t6, input_signal
    add $t6, $t6, $t5
    l.d $f4, 0($t6)

    mul.d $f6, $f2, $f4
    add.d $f0, $f0, $f6

    addi $sp, $sp, -8
    sw $t0, 0($sp)
    sw $t1, 4($sp)

    mov.d $f12, $f0
    jal round_double_1digit
    mov.d $f0, $f12

    lw $t0, 0($sp)
    lw $t1, 4($sp)
    addi $sp, $sp, 8

    addi $t1, $t1, 1
    j output_h_loop

save_output_y:
    sll $t3, $t0, 3
    la $t4, output_signal
    add $t4, $t4, $t3
    s.d $f0, 0($t4)

    addi $t0, $t0, 1
    j output_y_loop

output_done:
    lw $ra, 0($sp)
    addi $sp, $sp, 4
    jr $ra

# ----------------------------------------------------------------------------------
# -------------------------< Calculate MMSE Value >---------------------------------
# ----------------------------------------------------------------------------------

calculate_MMSE:
    addi $sp, $sp, -4
    sw $ra, 0($sp)

    li $t0, 0
    l.d $f0, zero_double

mmse_loop:
    bge $t0, 10, mmse_done

    sll $t1, $t0, 3

    la $t2, desired_signal
    add $t2, $t2, $t1
    l.d $f2, 0($t2)

    la $t3, output_signal
    add $t3, $t3, $t1
    l.d $f4, 0($t3)

    sub.d $f6, $f4, $f2

    addi $sp, $sp, -4
    sw $t0, 0($sp)

    mov.d $f12, $f6
    jal round_double_1digit
    mov.d $f6, $f12

    lw $t0, 0($sp)
    addi $sp, $sp, 4

    mul.d $f8, $f6, $f6

    addi $sp, $sp, -4
    sw $t0, 0($sp)

    mov.d $f12, $f8
    jal round_double_1digit
    mov.d $f8, $f12

    lw $t0, 0($sp)
    addi $sp, $sp, 4

    add.d $f0, $f0, $f8

    addi $sp, $sp, -4
    sw $t0, 0($sp)

    mov.d $f12, $f0
    jal round_double_1digit
    mov.d $f0, $f12

    lw $t0, 0($sp)
    addi $sp, $sp, 4

    addi $t0, $t0, 1
    j mmse_loop

mmse_done:
    l.d $f8, N_double
    div.d $f0, $f0, $f8

    addi $sp, $sp, -4
    sw $t0, 0($sp)

    mov.d $f12, $f0
    jal round_double_1digit
    mov.d $f0, $f12

    lw $t0, 0($sp)
    addi $sp, $sp, 4

    lw $ra, 0($sp)
    addi $sp, $sp, 4
    jr $ra

# ----------------------------------------------------------------------------------
# --------------------------< Write String To File >--------------------------------
# ----------------------------------------------------------------------------------

write_string_file:
    addi $sp, $sp, -8
    sw $ra, 0($sp)
    sw $a0, 4($sp)

    move $t0, $a0
    li $t1, 0

count_string_loop:
    lb $t2, 0($t0)
    beqz $t2, count_string_done
    addi $t1, $t1, 1
    addi $t0, $t0, 1
    j count_string_loop

count_string_done:
    li $v0, 15
    move $a0, $s6
    lw $a1, 4($sp)
    move $a2, $t1
    syscall

    lw $ra, 0($sp)
    addi $sp, $sp, 8
    jr $ra

# ----------------------------------------------------------------------------------
# --------------------------< Write Char To File >----------------------------------
# ----------------------------------------------------------------------------------

write_char_file:
    addi $sp, $sp, -4
    sw $ra, 0($sp)

    la $t0, one_char
    sb $a0, 0($t0)

    li $v0, 15
    move $a0, $s6
    la $a1, one_char
    li $a2, 1
    syscall

    lw $ra, 0($sp)
    addi $sp, $sp, 4
    jr $ra

# ----------------------------------------------------------------------------------
# --------------------------< Write Integer To File >-------------------------------
# ----------------------------------------------------------------------------------

write_int_file:
    addi $sp, $sp, -4
    sw $ra, 0($sp)

    move $t0, $a0
    bnez $t0, write_int_nonzero

    li $a0, 48
    jal write_char_file
    j write_int_done

write_int_nonzero:
    la $t1, int_buffer
    addi $t1, $t1, 19
    li $t2, 0

write_int_digit_loop:
    beqz $t0, write_int_output

    li $t3, 10
    div $t0, $t3
    mflo $t0
    mfhi $t4

    addi $t4, $t4, 48
    sb $t4, 0($t1)

    addi $t1, $t1, -1
    addi $t2, $t2, 1
    j write_int_digit_loop

write_int_output:
    addi $t1, $t1, 1

write_int_output_loop:
    beqz $t2, write_int_done

    lb $a0, 0($t1)
    jal write_char_file

    addi $t1, $t1, 1
    addi $t2, $t2, -1
    j write_int_output_loop

write_int_done:
    lw $ra, 0($sp)
    addi $sp, $sp, 4
    jr $ra

# ----------------------------------------------------------------------------------
# -----------------------< Write Double Round 1 Digit To File >----------------------
# ----------------------------------------------------------------------------------

write_double_1digit_file:
    addi $sp, $sp, -4
    sw $ra, 0($sp)

    l.d $f20, ten_double
    mul.d $f12, $f12, $f20

    round.w.d $f22, $f12
    mfc1 $t8, $f22

    bgez $t8, write_positive_1digit

    li $a0, 45
    jal write_char_file

    sub $t8, $zero, $t8

write_positive_1digit:
    li $t9, 10
    div $t8, $t9
    mflo $t6
    mfhi $t7

    move $a0, $t6
    jal write_int_file

    li $a0, 46
    jal write_char_file

    move $a0, $t7
    jal write_int_file

    lw $ra, 0($sp)
    addi $sp, $sp, 4
    jr $ra

# ----------------------------------------------------------------------------------
# --------------------< Print Double Round 1 Digit To Console >----------------------
# ----------------------------------------------------------------------------------

print_double_1digit_console:
    addi $sp, $sp, -4
    sw $ra, 0($sp)

    l.d $f20, ten_double
    mul.d $f12, $f12, $f20

    round.w.d $f22, $f12
    mfc1 $t8, $f22

    bgez $t8, print_positive_1digit

    li $a0, 45
    jal print_char_console

    sub $t8, $zero, $t8

print_positive_1digit:
    li $t9, 10
    div $t8, $t9
    mflo $t6
    mfhi $t7

    move $a0, $t6
    jal print_int_console

    li $a0, 46
    jal print_char_console

    move $a0, $t7
    jal print_int_console

    lw $ra, 0($sp)
    addi $sp, $sp, 4
    jr $ra

print_char_console:
    li $v0, 11
    syscall
    jr $ra

print_int_console:
    li $v0, 1
    syscall
    jr $ra
