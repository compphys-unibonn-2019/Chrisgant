#Hamiltonian function 1D. Gives matrix element from one state to another, look at test site (Wolfram.com) to check

H1<-function(initial, final,t,U,periodic = TRUE) {
    #interaction (static) case
	if(all(initial == final)){
		return(sum(initial[1:(length(initial)/2)] == initial[(length(initial)/2+1):(length(initial))] & initial[1:(length(initial)/2)] == 1) *U)
	}
	#kinetic (hopping) case
	if(sum(initial != final) == 2){
		if(sum(initial[1:(length(initial)/2)] != final[1:(length(initial)/2)]) != 1){
			w <- which(initial != final)
			d <- abs(w[1] - w[2])
			if(d == 1 | (periodic & d == (length(initial)/2-1)) ){
				n = sum(initial[w[1]:w[2]]) - 1
				return(-t * (-1)^n)
			}
		}else{
			#print("N_up and N_down not const")
		}
	}
	return(0)
}

#############Test area#############################
#Compare with: https://demonstrations.wolfram.com/HubbardModelInteractiveCalculatorFor1DSystems/
state1 = c(0,0,1,1,0,0,0,0)
state2 = c(1,0,1,0,0,0,0,0)

#print(H1(state1,state2,2,3,periodic = FALSE))