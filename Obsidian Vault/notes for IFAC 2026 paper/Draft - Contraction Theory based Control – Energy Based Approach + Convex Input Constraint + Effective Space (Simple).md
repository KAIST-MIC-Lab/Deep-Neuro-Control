## 1 Introduction 


## 2 Preliminaries
- 이론에 의하여 조건수로 정상상태 오차 그려짐

## 3 Conventional LMI–Approach (CV-STEM)
- 그와 관련된 LMI formulation
> Remark (Limitations) 수식 위주
> 수식적 한계 (미분항) + 입력포화
## 4 Main Result
- Effective space에 투영
	- contraction constraint
		- condition number minimization
			- angle objective function -> limitation (scale-invariant)
	- (convex) input saturation constraint 
		- norm constraint for demonstration
- 조건수는 scale-invariant -> energy based로 최소한의 입력이 가해지도록
	- energy constraint
- infeasible proble -> slack varialbe
> Remark about infeasiblity of energy constraint with respect to the error magnitude 
- Final formulation 
## 5 Numerical Validations
- Lorenz system
- disturbance
- saturation (norm-cosntraint)

- For comparative
	- (C1) Proposed one
	- (C2) CV-STEM
	- (C3) just ud
- C3가 가장 안 좋음
- C2는 포화가 되면서 됨
	- 조건수 1로 나오지만 contraction condition 위반함
	- 또한 포화됨
- C1은 모두 만족함
## 6 Conclusion
