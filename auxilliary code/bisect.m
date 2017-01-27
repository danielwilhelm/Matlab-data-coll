function x = bisect(a, b, f, maxIter)		% bisection algorithm to find root of f in [a,b]
	if (f(a)*f(b)>=0)
		warning('[bisect] f(a)*f(b)>=0');
		return;
	end;
	if a<b xl = a; xu = b; else xl = b; xu = a; end;
	if (f(xl)<0) incr=true; else incr=false; end;
		
	it = 1;
	while (it<=maxIter)
		fx = f((xu+xl)/2);
		if fx==0
			xl = (xu+xl)/2;
			xu = (xu+xl)/2;
		else
			if (incr)					
				if (fx > 0) xu = (xu+xl)/2; else xl = (xu+xl)/2; end;
			else
				if (fx > 0) xl = (xu+xl)/2; else xu = (xu+xl)/2; end;
			end;
		end;
		
		% fprintf('iteration %2i: %2.2f, %2.2f\n', it, xl, xu);
		it = it + 1;
	end;
	x = xu;
	
end