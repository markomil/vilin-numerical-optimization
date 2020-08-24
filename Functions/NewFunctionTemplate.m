  function [ outVal, outGr, outHes ] = POWER(x0, VGH )
% ================================================================
% Template for adding new multidimensional function to Vilin.
% To add new function modify this file and save to 'Functions/MultiDimensional/'.
% To add starting point for new function check 'Util/StartingPointGenerator.m'.
% ================================================================

    n = length(x0);

	outVal = 0;
    
    outHes = zeros(n,n);
    if VGH(1) > 0

        for i=1:n

			outVal = outVal + (x0(i)*i)^2;
        end
  
    end

    if VGH(2) > 0

	  for i=1:n
		outGr = zeros(n,1);
        outGr(i) = outGr(i) + 2*x0(i)*i;
      end
    end   

    if VGH(3) > 0

	  for i=1:n

			outHes = outHes + 2*i;

       end
     else

        outHes = 0;

     end    
