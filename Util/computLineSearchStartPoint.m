function lsStartPnt = computLineSearchStartPoint(fCurr, fPrev, iterNum, grad, dir, lsInitPnt)
% compute line search starting point according to Nocedal idea        

    if iterNum == 1 
        lsStartPnt = lsInitPnt;
    else
        lsStartPnt = abs(2*(fCurr-fPrev)/(grad'*dir));
        lsStartPnt = min(1, 1.01*lsStartPnt);
    end;
end
