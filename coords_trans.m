 function [out_coords]=coords_trans(input_coords,type)
 switch type
    case 'MNI305_2_MNI152'
    trans = [0.9975   -0.0073    0.0176   -0.0429
    0.0146    1.0009   -0.0024    1.5496
   -0.0130   -0.0093    0.9971    1.1840];
    input_coords = [input_coords,1]';
    out_coords = [trans*input_coords]';
 end
 end