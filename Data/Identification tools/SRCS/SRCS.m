clc
clear all
pkg load statistics

#file to load
load MeasuredData/AngularPendulum_Dif_L1.txt;
matrix = AngularPendulum_Dif_L1;
clear AngularPendulum_Dif_L1;

#properties
h = 0.01
iterations = 100
threshold = -0.001
countGenomes = 32
countExtSelection = 16
variableOfInteresst = 1
onlyMonomes =0

#dXdt: The vector which should be identified
dXdt = matrix(1:end-1,3);

#matrix: The trajectorie without dXdt 
matrix = matrix(:,2);
#matrix(:,3) = sqrt(matrix(:,1).*matrix(:,1) + matrix(:,2) .* matrix(:,2));

function EA(dataMatrix, dXdt, h, iterations, threshold, countGenomes, countExtSelection, variableOfInteresst, onlyMonomes)
  ea = EvolutionaryAlgorithm(size(dataMatrix,2), variableOfInteresst, countGenomes, countExtSelection, dataMatrix, dXdt, h, onlyMonomes);
  ea = fitnessEvaluation(ea, dataMatrix, 0); 
  loop(ea, iterations, threshold, 1);
 
endfunction

EA(matrix, dXdt, h, iterations, threshold, countGenomes, countExtSelection, variableOfInteresst, onlyMonomes);

