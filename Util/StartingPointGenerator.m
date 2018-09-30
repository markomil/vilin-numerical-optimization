function [ startingPoints ] = StartingPointGenerator(initialDimension)
% Generates default starting points for functions
% New points should be added at the end of this list
startingPoints = containers.Map();
startingPoints('GenRosenbrock') = repArray(initialDimension, -1.2, 1);
startingPoints('GenWhiteHolst') = repArray(initialDimension, -1.2, 1);
startingPoints('GenPSC1') = repArray(initialDimension, 3, 0.1);
startingPoints('ExtPSC1') = repArray(initialDimension, 3, 0.1);
startingPoints('ExtMaratos') = repArray(initialDimension, 1.1, 0.1);
startingPoints('ExtWhiteHolst') = repArray(initialDimension, -1.2, 1);
startingPoints('ExtRosenbrock') = repArray(initialDimension, -1.2, 1);
startingPoints('CUBE') = repArray(initialDimension, -1.2, 1);
startingPoints('Diagonal1') = oneNumber(1/initialDimension, initialDimension);
startingPoints('Diagonal2') = expArray(initialDimension, -1);
startingPoints('Diagonal3') = oneNumber(1, initialDimension);
startingPoints('Diagonal4') = oneNumber(1, initialDimension);
startingPoints('Diagonal5') = oneNumber(1.1, initialDimension);
startingPoints('Diagonal6') = oneNumber(1, initialDimension);
startingPoints('Diagonal9') = oneNumber(1, initialDimension);
startingPoints('DIXON3DQ') = oneNumber(-1, initialDimension);
startingPoints('ExtDENSCHNB') = oneNumber(1, initialDimension);
startingPoints('Hager') = oneNumber(1, initialDimension);
startingPoints('HIMMELHNew') = oneNumber(1.5, initialDimension);
startingPoints('Raydan1') = oneNumber(1, initialDimension);
startingPoints('Raydan2') = oneNumber(1, initialDimension);
startingPoints('ExtTridiag2') = oneNumber(1, initialDimension);
startingPoints('PertQuad') = oneNumber(0.5, initialDimension);
startingPoints('ExtPen') = aToB(1, initialDimension);
startingPoints('ExtHimmelblau') = oneNumber(1, initialDimension);
startingPoints('ExtTridiag1') = oneNumber(2, initialDimension);
startingPoints('GenTridiag1') = oneNumber(2, initialDimension);
startingPoints('ExtPowell') = repArray(initialDimension, 3,-1,0,1);
startingPoints('FullHessian2') = oneNumber(0.01, initialDimension);
startingPoints('ExtBD1') = oneNumber(0.1, initialDimension);
startingPoints('QuadQF1') = oneNumber(1, initialDimension);
startingPoints('QUARTC') = oneNumber(2, initialDimension);
startingPoints('ExtQuadPenQP1') = oneNumber(1, initialDimension);
startingPoints('ExtHiebert') = oneNumber(0, initialDimension);
startingPoints('QuadQF2') = oneNumber(0.5, initialDimension);
startingPoints('FLETCHCR') = oneNumber(0, initialDimension);
startingPoints('LIARWHDNew') = oneNumber(4, initialDimension);
startingPoints('BDQRTIC') = oneNumber(1, initialDimension);
startingPoints('TRIDIA') = oneNumber(1, initialDimension);
startingPoints('ARGLINB') = oneNumber(1, initialDimension);
startingPoints('NONDIA') = oneNumber(-1, initialDimension);
startingPoints('NONDQUAR') = repArray(initialDimension, 1, -1);
startingPoints('NONSCOMP') = oneNumber(3, initialDimension);
startingPoints('DQDRTIC') = oneNumber(3, initialDimension);
startingPoints('EG2') = oneNumber(1, initialDimension);
startingPoints('SINENew') = oneNumber(1, initialDimension);
startingPoints('SINQUADNew') = oneNumber(0.1, initialDimension);
startingPoints('PartPertQuad') = oneNumber(0.5, initialDimension);
startingPoints('AlmostPertQuad') = oneNumber(0.5, initialDimension);
startingPoints('PertTridiagQuad') = oneNumber(0.5, initialDimension);
startingPoints('Staircase1New') = oneNumber(1, initialDimension);
startingPoints('Staircase2New') = oneNumber(2, initialDimension);
startingPoints('New_ExtFreudAndRoth') = repArray(initialDimension, 0.5, -2);
startingPoints('New_ExtBeale') = repArray(initialDimension, 1, 0.8);
startingPoints('New_ExtTET') = oneNumber(2,initialDimension);
startingPoints('New_FullHessian1') = oneNumber(0.01,initialDimension);
startingPoints('New_Staircase1') = oneNumber(1,initialDimension);

% starting point for custom function
% startingPoints('CustomFunction') = point_generator(initialDimension...);

% point_generators:

function sp = oneNumber(number, dimension)
% returns number repeted dimesion times
sp = repmat(number, 1, dimension);

function sp = aToB(a, b)
% return integers between a and b
sp = linspace(a, b, b - a + 1);

function sp = expArray(dimension, exp)
% returns n^exp array with dimension elements
sp = zeros(1, dimension);
for i=1:dimension
    sp(i) = round(i^exp, 5, 'significant');
end

function sp = repArray(dimension, varargin)
% returns array with repeted provided elements of lentght dimension
% if no elements are provided returns zeros(1, dimension)
if nargin > 0
    aux = cell2mat(varargin);
    sp = repmat(aux, 1, round(dimension / length(aux)) + 1);
    sp = sp(1:dimension);
else
    sp = zeros(1, dimension);
end
