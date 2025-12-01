  - M이 있음
	- 수축조건이 BMI로 나옴; M기준
	- 이걸 특정 벡터(실제어에 사용하는 e)에 투영해서 쉽게 만듬

  
- M을 바로 푸는것 보다 M =: m*I + \Delta M으로 분리
    - 목적함수를 min cond(M)임으로 \Dleta_M이 줄어들도록
    - m은 제어오차와 수축률에 비례하도록 되면 좋겠음
- 목적함수 정의
    - m기준으로
        - 스칼라문제 -> 뒤집어진 이차함수 -> 바이어스(미분항) 무시 -> 0을 지나는 이차함수
        - 그중 하나만 고르면 되는데 m만 있다면 어디든 조건수가 1인 지점 (이차함수가 음인 지점에서)
        - 이중 가장 적은 제어입력인 m을 사용하면 좋음
            - 목적함수를 이차함수의 0이 아닌 근에 대한 이차식으로 작성

        - ->>> 조건수가 1인 M이 있다고 가정하며 상수배라고 할때

    - \Delta M포함
        - 목적함수를 위와 비슷한 접근 사용
            - 바이어스 무시 -> 0을 지나는 타원
            - 조건수가 1인 최적점 수축조건이 0이되는 지점
                - 수축조건이 0이 되도록 작성
            - 이에 더불어서 \Delta M이 0이 되도록

        - ->>> 실제 M 기준으로 확장함

    - 포화 조건
        - 알아서 포함


- [NEW] Contraction based control
    - better calculation of contraction metric [IFAC]
    - input saturation [IFAC] -> DL
    - (error) state saturation [CDC Feb]
    - ->> jour.
- FUZZY -> DL (offline, slow, global) + NAC (online, fast, local)
    - ->> conf,
    - DL: contraction, data-driven
    - NAC: CONAC