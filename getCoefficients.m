% generate sequence of regression coefficients
function coeff = getCoefficients(M, model, scale)

	switch model

		case 'exp'
			coefftmp = exp(-linspace(1,M,M));
			coeff = 0.5*scale*coefftmp'*10;
		
		case 'power'
			coefftmp = 2.^(-linspace(1,M,M));
			coeff = 0.5*scale*coefftmp';

		case 'geom'
			coefftmp = linspace(1,M,M).^(-3);
			coeff = 0.5*scale*coefftmp';

		case 'sparse'
			coefftmp = [ones(1,3) zeros(1,M-3)];
			coeff = 0.5*scale*coefftmp';

		case 'lin-exp'
			coeff = 0.5*scale*[ linspace(3,1,5) exp(-linspace(1,M-5,M-5)) ]';
		
		case 'lin-power'
			coeff = 0.5*scale*[ linspace(3,1,5) 2.^(-linspace(1,M-5,M-5)) ]';

		case 'lin-geom'
			coeff = 0.5*scale*[ linspace(3,1,5) linspace(2,M-5,M-5).^(-3) ]';	

		case 'lin-sparse'
			coeff = 0.5*scale*[linspace(3,1,5) zeros(1,M-5)]';
	end

end