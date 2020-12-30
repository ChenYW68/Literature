rm(list=ls())

require("grplasso")
require("splines")
require("base")
require("psych")
require("MASS")
require("splus2R")
require("mnormt")
require("nlme")
require("xtable")
require("matrixcalc")


Design.M = function(EX,J,m,Seed){

	# Arguments -
	# EX is the for two different kind of examples of simulating the covariates
	# J is total number of additive covariates in the model
	# m is the lattice size of the square lattice
	# Seed controls the simulation seed based on which MC simulation we are in
	
	# Returns - 
	# design matrix as "des" and indicator for true beta as "itb"
	
	n = m^2			
	# n is the Sample size

	if (EX == 1){

		des = matrix(rep(0,n*J),byrow=TRUE,ncol=J)
		for(j in 1:J){

			Seed = Seed + 231
			set.seed(Seed)
			des[,j] = runif(n,0,1)

		}

		f1 = 3*(2*des[,2] - 1)^2
		f2 = (2*sin(2*pi*des[,3]))/(2 - sin(2*pi*des[,3]))

		q = 2
		des = cbind(f1,f2,des[,(q+1):J])

		names1 = paste("f",1:q)
		names2 = paste("X",(q+1):J)
		colnames(des) = c(names1,names2)

		itb = c(rep(1,q),rep(0,J-q))

	}
	
	if (EX == 2){

		W = matrix(rep(0,n*J),byrow=TRUE,ncol=J)
		for(j in 1:J){

			Seed  = Seed + 231
			set.seed(Seed)
			W[,j] = runif(n,0,1)

		}

		Seed = Seed + 457
		set.seed(Seed)
		U = runif(n,0,1)
		t = 3

		des = matrix(rep(0,n*J),byrow = TRUE,ncol = J)
		for(j in 1:J){

			des[,j] = (W[,j] + t*U)/(1 + t)

		}

		f1 = 5*des[,1]
		f2 = 3*(2*des[,2] - 1)^2
		f3 = (2*sin(2*pi*des[,3]))/(2 - sin(2*pi*des[,3]))
		f4 = 0.6*sin(2*pi*des[,4]) + 1.2*cos(2*pi*des[,4])+1.8*(sin(2*pi*des[,4]))^2 
			+ 2.4*(cos(2*pi*des[,4]))^3 + 3.0*(sin(2*pi*des[,4]))^3

		q = 4
		des = cbind(f1,f2,f3,f4,des[,(q+1):J])

		names1 = paste("f",1:q)
		names2 = paste("X",(q+1):J)
		colnames(des) = c(names1,names2)

		itb = c(rep(1,q),rep(0,J-q))

	}

	return(list(des=des,itb=itb))

}



Cov.M_Spatial.W = function(par,type,m){ 

	# Arguments - 
	# par is the value of the parameter in the covariance model
	# type is the type of covariance model
	# m is the lattice size

	# Exp = Exponential
	# Mat_3_2 = Matern_3/2
	# Mat_5_2 = Matern_5/2
	# Gauss = Gaussian
	# InvMQ = Inverse MQ

	# Returns - 
	# both covariance matrix and the spatial weight matrix

	n = m^2

	M = matrix(1:m,m,m)
	Ind = cbind(c(t(M)),c(M))

	G = matrix(rep(1,n^2),byrow=TRUE,ncol=n)
	
	for (i in 1:(n-1)){

		for (j in (i+1):n){
	
			tmp = vecnorm(Ind[i,]-Ind[j,])
			if (type == "Exp") { G[i,j] = exp(-par*tmp) }
			else
			if (type == "Mat_3_2") { G[i,j] = (1 + (sqrt(3)*tmp)/par)*
								exp(-(sqrt(3)*tmp)/par) }
			else
			if (type == "Mat_5_2") { G[i,j] = (1 + (sqrt(5)*tmp)/par + 
								((sqrt(5)*tmp)/par )^2/3)*
								exp(-(sqrt(5)*tmp)/par) }
			else
			if (type == "Gauss") { G[i,j] = exp(-(par^2)*(tmp^2)) }
			else
			if (type == "InvMQ") { G[i,j] = (1 + tmp^2 )^(-par)} 
			else { G[i,j] = 0 }

			G[j,i] = G[i,j]

		}

	}

	return(G)

}



M.name = c('Exp','Mat_3_2','Mat_5_2','Gauss')
W.name = c('I','Gauss','InvMQ','True')

PAR = array(0,c(6,2))
PAR[1,] = c(0.5,1)	# for Exp = Exponential
PAR[2,] = c(2.5,1.5)	# for Mat_3/2 = Matern_3/2
PAR[3,] = c(2.5,1.5)	# for Mat_5/2 = Matern_5/2
PAR[4,] = c(1.5,2.5)	# for Gauss = Gaussian.model
PAR[5,] = c(2,3.5)	# for Gauss = Gaussian.weight
PAR[6,] = c(1.5,2.5)	# for InvMQ = Inverse MQ


M.par = matrix(rep(0,2*length(M.name)),byrow = TRUE,ncol = 2)
M.par[1,] = PAR[1,]
M.par[2,] = PAR[2,]
M.par[3,] = PAR[3,]
M.par[4,] = PAR[4,]

W.par = matrix(rep(0,2*length(W.name)),byrow = TRUE,ncol = 2)
W.par[1,] = c(0,0)
W.par[2,] = PAR[5,]
W.par[3,] = PAR[6,]
W.par[4,] = c(0,0)

Gauss.par.seq = seq(0.5,5,0.25)
InvMQ.par.seq = seq(0.25,5,0.25)


ave_sd.sel.var = function(MC,m,J,EX,Model,Weight,par.m,par.w,seednum) {

	# Arguments - 
	# MC index of Monte Carlo simulation
	# m is the size of the lattice, n = m*m is the sample size
	# J is the number of total parameters
	# EX the example we will be using, the two setups in out paper
	# Model is the true covariance model using which the data wil be generated
	# Weight is the spatial weight matrix used for weighed objective function
	# seednum is used for the index to change design matrix of different simulation sets

	# Returns -
	# Monte carlo mean and standard deviation from several simulations
	# for Group LASSO and adaptive group LASSO and also True and False Positive indices 
	
	# Note on penalty parameter -
	# The optimal penalty parameter selection is cross validated over a range of parameters from spatial weight matrix

	GL0 = rep(0,MC)		# Vector to store mean non-zero variables selected through GLASSO
	AGL0 = rep(0,MC)		# Vector to store mean non-zero variables selected through AGLASSO
	TP.GL = rep(0,MC)		# Vector to store mean True Positive through GLASSO
	FP.GL = rep(0,MC)		# Vector to store mean False Positive through GLASSO
	TP.AGL = rep(0,MC)	# Vector to store mean True Positive through AGLASSO
	FP.AGL = rep(0,MC)	# Vector to store mean False Positive through AGLASSO

	for(iS in 1:MC) {

		n = m^2
		seednum = seednum + 111
		set.seed(seednum)

		Dsgn = Design.M(EX,J,m,iS+seednum) 
		snr = 3							## Signal to noise ratio
		D.M = Dsgn$des						## Temporary variable for Design matrix
		sd_f = sqrt(tr(var(D.M)))				## Measuring Signal for additive models
		sigma.sq.true =(sd_f/snr)^2				## True value of sigma square
		Sigma.true = sigma.sq.true*Cov.M_Spatial.W(par.m,Model,m)

		seednum = 2*(seednum + 111)
		set.seed(seednum)

		E = t(rmvnorm(1,rep(0,n),Sigma.true))
		z = D.M
		bs.x = NULL
		M.n = rep(0,J)

		for(j in 1:J){

			Q0 = quantile(z[,j],c(0.25,0.5,0.75))	## Interior knots
			temp1 = bs(z[,j],knots = Q0,Boundary.knots = range(z[,j]))
			temp2 = ncol(temp1)
			bs.x = cbind(bs.x,temp1)
			M.n[j] = temp2

		}

		mu = rep(0,n) # WLOG we are taking mu = 0

		y = mu + D.M[,1] + D.M[,2] + D.M[,3] + D.M[,4] + E 
		par.w = seq(0.5,5,0.25)
		lambda = rep(0,length(par.w))
		beta.grpLasso = matrix(rep(0, sum(M.n)*length(par.w) ), byrow = T, ncol = length(par.w) )
		EGIC = rep(0,length(par.w))

		for (l in 1:length(par.w)) {

			Sigma.work = sigma.sq.true*Cov.M_Spatial.W(par.w[l],Weight,m)
			tmp0 = chol(Sigma.work)
			RW = t(solve(tmp0))
			spnormRW = sqrt(1/min(eigen(Sigma.work)$values))
	
			alpha = 0 
			X = RW%*%(bs.x)
			Y = RW%*%y
			C = 2.5
			lambdamin1 = C*spnormRW*sqrt((n^(1+alpha))*M.n[1]*log(J*M.n[1]))
			lambda[l] = lambdamin1

			index = c(rep(seq(1:length(M.n)),M.n)) 		# Grouping of the variables
	
			fit = grplasso(X,Y,index = index,lambda = lambdamin1,coef.init = c(rep(0,dim(X)[2])),
				model = LinReg(),center = FALSE,standardize = TRUE,control = grpl.control(update.hess = "lambda",
				trace = 0))

			beta.g = as.matrix(fit$coeff)
			beta.grpLasso[,l] = beta.g

			c.M = c(0,cumsum(M.n))
			eta = rep(0,J)

			for(j in 1:J){

				eta[j] = vecnorm(beta.g[(c.M[j]+1):c.M[j+1]])

			}

			df.lambda = sum(sign(eta))*M.n[1]
			EGIC[l] = log (  vecnorm(Y - predict(fit))^2 ) + df.lambda*(log(log(n))/n)

		}

		lambda.opt = lambda[which.min(EGIC)]
		beta.grpLasso.opt = beta.grpLasso[,which.min(EGIC)]

		c.M = c(0,cumsum(M.n))
		beta.gL = matrix(rep(0,J*M.n[1]),byrow = TRUE,ncol = J)
		eta = rep(0,J)

		for(j in 1:J){

			beta.gL[,j] = beta.grpLasso.opt[(c.M[j]+1):c.M[j+1]]
			eta[j] = vecnorm(beta.gL[,j])

		}

		colnames(beta.gL) = c(paste("Group",1:J))
		ieb.gl = sign(eta)
		pos.1 = which(ieb.gl == 1)

		itb = Dsgn$itb
		itb.tem = sign(abs(itb-1))

		TP.GL[iS] = sum(itb*ieb.gl)
		FP.GL[iS] = sum(itb.tem*ieb.gl)

		##################################################################
		##################################################################
		###										   ###
		### To find the set of components whose group LASSO L2 norm is ###
		###                        greater than 0                      ###
		###										   ###
		##################################################################
		##################################################################

		A0.tilde = NULL
	
		for(j in 1:J){

			if (eta[j] == 0) { A0.tilde = c(A0.tilde,j) } else  { A0.tilde = A0.tilde }

		}


		A = 1:J						## Index for set of all the components
		A1.tilde = A[is.na(pmatch(A,A0.tilde))]	
		## Index set for components whose L2 norm is greater than zero

		GL0[iS] = length(A1.tilde)

		if (GL0[iS] > 0){

			#############################################
			#############################################
			##							 ##
			##           Adaptive group LASSO          ##
			##							 ##
			#############################################
			#############################################

			z1 = as.matrix(D.M[,A1.tilde])
			eta.f = eta[A1.tilde]
			bs.x.new = NULL
			M.n1 = rep(0,GL0[iS])

			for(j in 1:GL0[iS]){

				Q1 = quantile(z1[,j],c(0.25,0.5,0.75))	## Interior knots
				temp3 = bs(z1[,j],knots = Q1,Boundary.knots = range(z1[,j]))/sqrt(eta.f[j])
				temp4 = ncol(temp3)
				bs.x.new = cbind(bs.x.new,temp3)
				M.n1[j] = temp4

			}

			lambda2 = rep(0,length(par.w))
			EGIC.agl = rep(0,length(par.w))
			beta.AgrpLasso = matrix(rep(0, sum(M.n1)*length(par.w) ), byrow = T, ncol = length(par.w) )

			for (l in 1:length(par.w)) {

				X.new = RW%*%(bs.x.new)
				lambdamin2 = min( (sqrt(M.n1[1])*lambdamin1) , (n/M.n1[1]) )
				lambda2[l] = lambdamin2

				index.agl = c(rep(seq(1:length(M.n1)),M.n1)) # Grouping of the variables

				fit.agl = grplasso(X.new,Y,index = index.agl,lambda = lambdamin2,coef.init = c(rep(0,dim(X.new)[2])),
					model = LinReg(),center = FALSE,standardize = TRUE,control = grpl.control(update.hess = "lambda",
					trace = 0))

				beta.Ag = as.matrix(fit.agl$coeff)
				beta.AgrpLasso[,l] = beta.Ag

				c.M1 = c(0,cumsum(M.n1))
				eta.AgL = rep(0,GL0[iS])

				for(j in 1:GL0[iS]){

					eta.AgL[j] = vecnorm(beta.Ag[(c.M1[j]+1):c.M1[j+1]])

				}

				df.lambda.agl = sum(sign(eta.AgL))*M.n1[1]
				EGIC.agl[l] = log (  vecnorm(Y - predict(fit.agl))^2 ) + df.lambda.agl*(log(log(n))/n)

			}
		
			lambda.agl.opt = lambda2[which.min(EGIC.agl)]
			beta.AgrpLasso.opt = beta.AgrpLasso[,which.min(EGIC.agl)]

			c.M1 = c(0,cumsum(M.n1))
			beta.AgL = matrix(rep(0,GL0[iS]*M.n1[1]),byrow = TRUE,ncol = GL0[iS])
			eta.AgL = rep(0,GL0[iS])

			for(j in 1:GL0[iS]){

				beta.AgL[,j] = beta.AgrpLasso.opt[(c.M1[j]+1):c.M1[j+1]]
				eta.AgL[j] = vecnorm(beta.AgL[,j])

			}

			colnames(beta.AgL) = c(paste("Group",1:GL0[iS]))

			ieb.agl0 = sign(eta.AgL)			
			ieb.agl = rep(0,J)
			ieb.agl[pos.1] = ieb.agl0

			TP.AGL[iS] = sum(itb*ieb.agl)
			FP.AGL[iS] = sum(itb.tem*ieb.agl)

			##################################################################
			##################################################################
			###										   ###
			### To find the set of components whose group LASSO L2 norm is ###
			###                        greater than 0                      ###
			###										   ###
			##################################################################
			##################################################################

			A2.tilde = NULL

			for(j in 1:GL0[iS]){

				if (eta.AgL[j] == 0) { A2.tilde = c(A2.tilde,j) } else  { A2.tilde = A2.tilde }

			}


			A = 1:GL0[iS]		## Index for set of all the components
			A3.tilde = A[is.na(pmatch(A,A2.tilde))] 
			## Index set for components whose L2 norm is greater than zero

			AGL0[iS] = length(A3.tilde)

		}
		
		if (GL0[iS] == 0){

			lambdamin2 = 'N.A'
			AGL0[iS] = 0

		}
		

	}

	ret1 = mean(GL0)
	ret2 = sd(GL0)
	
	ret3 = mean(TP.GL)
	ret4 = sd(TP.GL)
	
	ret5 = mean(FP.GL)
	ret6 = sd(FP.GL)

	ret7 = mean(AGL0)
	ret8 = sd(AGL0)

	ret9 = mean(TP.AGL)
	ret10 = sd(TP.AGL)
	
	ret11 = mean(FP.AGL)
	ret12 = sd(FP.AGL)

	return(list(gL.mean = ret1,gL.sd = ret2,true.pos_gL.mean = ret3,true.pos_gL.sd = ret4,false.pos_gL.mean = ret5,false.pos_gL.sd = ret6,lambdamin1 = lambda.opt,
			agL.mean = ret7,agL.sd = ret8,true.pos_agL.mean = ret9,true.pos_agL.sd = ret10,false.pos_agL.mean = ret11,false.pos_agL.sd = ret12,lambdamin2 = lambda.agl.opt))

}


ave_sd.sel.var(MC = 2,m = 12,J = 15,EX = 2,Model = 'Exp',Weight = 'Gauss',par.m = M.par[1,2],par.w = seq(0.5,5,0.25),seednum = 101)


