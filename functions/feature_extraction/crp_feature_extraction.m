function [ feature_struct ] = crp_feature_extraction( cell_data, h )

    [ feature_S_values ] = crpSummary( h,cell_data );
    %    feature_mu_C = z1([on_C_mu off_C_mu gray_C_mu crp_C_mu amp_C_mu]);
    %     feature_vec_C = [feature_mu_C;[ on_C_std  off_C_std  gray_C_std  crp_C_std  amp_C_std]];
        % constants
    crp_part = [13    20    26    60]-2; % 5 states extracted with crp_partition_scrp
    amp_part = [16    27    35    48    57    71    79    90   101   112   122]-2; % 11 states extracted with crp_partition_scrp
    contrast_space = linspace(0,255,133);
    freq_space = linspace(0.75,6,61);
    freqs = freq_space(crp_part);
    contrasts = round(contrast_space(amp_part));
    
    feature_struct.S_on_mu = feature_S_values(1);
    feature_struct.S_on_std = feature_S_values(2);
    feature_struct.S_off_mu = feature_S_values(3);
    feature_struct.S_off_std = feature_S_values(4);
    feature_struct.S_gray_mu = feature_S_values(5);
    feature_struct.S_gray_std = feature_S_values(6);
    
    feature_struct.S_crp_mu_075 = feature_S_values(7);
    feature_struct.S_crp_std_075 = feature_S_values(8);
    feature_struct.S_crp_mu_162 = feature_S_values(9);
    feature_struct.S_crp_std_162 = feature_S_values(10);
    feature_struct.S_crp_mu_223 = feature_S_values(11);
    feature_struct.S_crp_std_223 = feature_S_values(12);
    feature_struct.S_crp_mu_276 = feature_S_values(13);
    feature_struct.S_crp_std_276 = feature_S_values(14);
    feature_struct.S_crp_mu_573 = feature_S_values(15);
    feature_struct.S_crp_std_573 = feature_S_values(16);
    
    feature_struct.S_amp_mu_0 = feature_S_values(17);
    feature_struct.S_amp_std_0 = feature_S_values(18);
    feature_struct.S_amp_mu_25 = feature_S_values(19);
    feature_struct.S_amp_std_25 = feature_S_values(20);
    feature_struct.S_amp_mu_46 = feature_S_values(21);
    feature_struct.S_amp_std_46 = feature_S_values(22);
    feature_struct.S_amp_mu_62 = feature_S_values(23);
    feature_struct.S_amp_std_62 = feature_S_values(24);
    feature_struct.S_amp_mu_87 = feature_S_values(25);
    feature_struct.S_amp_std_87 = feature_S_values(26);
    feature_struct.S_amp_mu_104 = feature_S_values(27);
    feature_struct.S_amp_std_104 = feature_S_values(28);
    feature_struct.S_amp_mu_131 = feature_S_values(29);
    feature_struct.S_amp_std_131 = feature_S_values(30);
    feature_struct.S_amp_mu_147 = feature_S_values(31);
    feature_struct.S_amp_std_147 = feature_S_values(32);
    feature_struct.S_amp_mu_168 = feature_S_values(33);
    feature_struct.S_amp_std_168 = feature_S_values(34);
    feature_struct.S_amp_mu_189 = feature_S_values(35);
    feature_struct.S_amp_std_189 = feature_S_values(36);
    feature_struct.S_amp_mu_211 = feature_S_values(37);
    feature_struct.S_amp_std_211 = feature_S_values(38);
    feature_struct.S_amp_mu_230 = feature_S_values(39);
    feature_struct.S_amp_std_230 = feature_S_values(40);
end