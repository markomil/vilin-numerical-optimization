function [ fmin, xmin, iterNum, cpuTime, evalNumbers, valuesPerIter ] = constKorak( f, objekat )
%Constant step scanning

 tic;

maxIter = objekat.max_iteration_no;
evalNumbers = EvaluationNumbers(0,0,0);
valuesPerIter = PerIteration(maxIter);
 a = objekat.od;
 b = objekat.do;
 delta = objekat.step_size;

fmin=feval(f,a);
evalNumbers.incrementBy([1 0 0]);
iterNum = 1;
valuesPerIter.setFunctionVal(iterNum, fmin);
xmin=a;
x=a;

while x<=b
    iterNum = iterNum + 1;
    x=x+delta;
    ftr=feval(f,x);
    evalNumbers.incrementBy([1 0 0]);
    valuesPerIter.setFunctionVal(iterNum, ftr);
    if ftr<fmin
        fmin=ftr;
        xmin=x;
    end
end

cpuTime = toc;
