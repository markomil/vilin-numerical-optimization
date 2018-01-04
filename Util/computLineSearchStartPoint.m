function lsStartPnt = computLineSearchStartPoint(fCurr, fPrev, grad, dir)
% compute line search starting point according to the idea proposed by
% Nocedal and Wright 'Numerical optimization'
% lecture: 'The initial step length', page 58, chapter 3

        lsStartPnt = abs(2*(fCurr-fPrev)/(dir*grad));
        lsStartPnt = min(1, 1.01*lsStartPnt);
    
end
