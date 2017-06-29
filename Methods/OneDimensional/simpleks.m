function [ fmin, xmin, iterNum, cpuTime, evalNumbers, valuesPerIter ] = simpleks( f, objekat )
%Onedimensional simplex method

tic;

maxIter = objekat.max_iteration_no;
evalNumbers = EvaluationNumbers(0,0,0);
valuesPerIter = PerIteration(maxIter);
a = objekat.starting_point;
delta = objekat.step_size;
deltaMin = objekat.step_size_min;

if objekat.variables_no>length(a)
    a=repmat(a, 1, objekat.variables_no);
end

fmin=feval(f,a);
evalNumbers.incrementBy([1 0 0]);
xmin=a;
x = a;
L = 4;
fpr=fmin;
iterNum=1;

while abs(delta) > deltaMin && iterNum < maxIter
    x=x+delta;
    ftr=feval(f,x);
    evalNumbers.incrementBy([1 0 0]);
    valuesPerIter.setFunctionVal(iterNum, ftr);

    if ftr<=fmin
        fmin=ftr;
        xmin=x;
    end

    if ftr>fpr
        delta=-delta/L;
    end

    iterNum=iterNum+1;
    fpr=ftr;
end

cpuTime=toc;

valuesPerIter.trim(iterNum);

