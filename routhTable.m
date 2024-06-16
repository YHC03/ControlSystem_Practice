%routhTable(polesArray) Function
%
% Creates Routh Table and print Stability of the Function based on the Calculation of Routh Table
%
% polesArray: array for poles
%
% Written By: YHC03
function outputArray = routhTable(polesArray)

	%evenLocale: 모든 행이 0인 Row의 위치
	evenLocale = 0;
	% routhValue: Routh Table to Output
	routhValue = [];
	% routhRow: A Row for routhValue. For Creation of a Row.
	routhRow = [];

	% Pole의 차수 구하기
	polesLength = length(polesArray) - 1;
	% 행렬의 Column 길이 구하기
	polesLengthEven = double(int32(polesLength / 2) + 1 - rem(polesLength, 2));


	% 행렬의 첫째 Row의 값 구하기
	for i = 1:polesLengthEven
		routhRow = [routhRow, polesArray(2 * i - 1)];
	end
	% routhValue에 결과 반영
	routhValue = [routhValue; routhRow];


	% routhRow 초기화
	routhRow = [];
	% 행렬의 두번째 Row의 값 구하기
	for i = 1:(polesLengthEven - 1 + rem(polesLength, 2))
		routhRow = [routhRow, polesArray(2 * i)];
	end
	% 입력된 행렬의 항의 갯수가 홀수인 경우, 행렬 마지막에 0 추가
	if rem(polesLength, 2) == 0
		routhRow = [routhRow, 0];
	end
	
	% 행렬의 첫번째 값이 0인 경우를 찾는다
	if routhRow(1) == 0
		for j = 2:length(routhRow)
			if routhRow(j) ~= 0
				% 행렬의 첫번째 값만 0인 경우
				fprintf("Found 0 at leftmost side of Calculation\n");
				disp([routhValue; routhRow]);
				fprintf("...\n");
				outputArray = [];
				return;
			end
		end
		
		% 전체 행렬이 0인 경우
		fprintf("Row of Zero was found at s^%d\n", polesLength - 1);
		evenLocale = 1; % Row of Zero의 위치 기록
		routhRow = []; % routhRow 초기화

		% 직전 값을 미분하여 routhRow 행렬을 채운다
		for targetRouthColumn = 1:polesLengthEven
			if (polesLength - (targetRouthColumn - 1) * 2) > 0
				% 해당 행렬의 위치의 차수가 양수인 경우
				routhRow = [routhRow, routhValue(1, targetRouthColumn) * (polesLength - (targetRouthColumn - 1) * 2)];
			else
				% 해당 행렬의 위치의 차수가 0이나 음수인 경우
				routhRow = [routhRow, 0];
			end
		end
	end

	% 결과 행렬에 해당 값 추가
	routhValue = [routhValue; routhRow];


	% 나머지 Row에 대한 연산 수행
	for prevYCurr = 2:polesLength
		routhRow = []; % routhRow 초기화

		% 해당 Row에 대한 연산을 수행한다
		for xCurr = 2:(polesLengthEven + 1)
			if xCurr >= (polesLengthEven + 1)
				% 행렬 외부를 참조하려 하는 경우, 계산값은 무조건 0이다.
				routhRow = [routhRow, 0];
			else
				% Routh Table 연산을 수행한다
				routhRow = [routhRow, -1 * (routhValue(prevYCurr - 1, 1) * routhValue(prevYCurr, xCurr) - routhValue(prevYCurr - 1, xCurr) * routhValue(prevYCurr, 1)) / routhValue(prevYCurr, 1)];
			end
		end
	
		% 행렬의 첫번째 값이 0인 경우를 찾는다
		if routhRow(1) == 0
			for j = 2:length(routhRow)
				if routhRow(j) ~= 0
					% 행렬의 첫번째 값만 0인 경우
					fprintf("Found 0 at leftmost side of Calculation\n");
					disp([routhValue; routhRow]);
					fprintf("...\n");
					outputArray = [];
					return;
				end
			end
		
			% 전체 행렬이 0인 경우
			fprintf("Row of Zero was found at s^%d\n", polesLength - 1);
			evenLocale = prevYCurr; % Row of Zero의 위치 기록
			routhRow = []; % routhRow 초기화

			% 직전 값을 미분하여 routhRow 행렬을 채운다
			for targetRouthColumn = 1:polesLengthEven
				if (polesLength - prevYCurr + 1 - (targetRouthColumn - 1) * 2) > 0
					% 해당 행렬의 위치의 차수가 양수인 경우
					routhRow = [routhRow, routhValue(prevYCurr, targetRouthColumn) * (polesLength - prevYCurr + 1 - (targetRouthColumn - 1) * 2)];
				else
					% 해당 행렬의 위치의 차수가 0이나 음수인 경우
					routhRow = [routhRow, 0];
				end
			end
		end

		% 결과 행렬에 해당 값 추가
		routhValue = [routhValue; routhRow];
	end

	
	% Root 판별
	% 임시 변수
	prevSign = 0; % 직전에 읽은 값의 부호(1: 양수, -1: 음수. 0은 나타나지 않음)
	signChanges = 0; % 부호 변화 횟수

	% 결괏값 저장 변수
	non_JW_at_P = 0; % Even에서 jw-Axis에 존재하지 않는 해의 갯수
	jw_at_P = 0; % Even에서 jw-Axis에 존재하는 해의 갯수
	rhp_at_Q = 0; % Odd에서 RHP에 존재하는 해의 갯수
	lhp_at_Q = 0; % Odd에서 LHP에 존재하는 해의 갯수

	if evenLocale > 0 % Row of Zero가 나타난 경우
		% 최초의 부호값 부여
		if routhValue(1, 1) > 0
			prevSign = 1;
		else
			prevSign = -1;
		end

		% 최초의 Row of Zero가 나타난 위치 직전까지 진행
		for i = 2:evenLocale
			if routhValue(i, 1) > 0 && prevSign == -1
				% 음수에서 양수로 변한 경우
				signChanges = signChanges + 1;
				prevSign = 1;
			elseif routhValue(i, 1) < 0 && prevSign == 1
				% 양수에서 음수로 변한 경우
				signChanges = signChanges + 1;
				prevSign = -1;
			end
		end

		% 부호의 변화 횟수는 Odd에서 RHP에 존재하는 해의 갯수
		% 위의 반복 횟수에서 나머지 횟수는 Odd에서 LHP에 존재하는 해의 갯수
		rhp_at_Q = signChanges;
		lhp_at_Q = evenLocale - 1 - signChanges;


		% 나머지 구간 연산
		signChanges = 0;
		% 구간 기준, 최초의 부호값 부여
		if routhValue(evenLocale + 1, 1) > 0
			prevSign = 1;
		else
			prevSign = -1;
		end

		% 최초의 Row of Zero가 나타난 위치 직후부터 끝까지 진행
		for i = (evenLocale + 2):(polesLength + 1)
			if routhValue(i, 1) > 0 && prevSign == -1
				% 음수에서 양수로 변한 경우
				signChanges = signChanges + 1;
				prevSign = 1;
			elseif routhValue(i, 1) < 0 && prevSign == 1
				% 양수에서 음수로 변한 경우
				signChanges = signChanges + 1;
				prevSign = -1;
			end
		end

		% 부호의 변화 횟수는 Even에서 RHP와 LHP에 존재하는 해의 갯수
		% 위의 반복 횟수에서 나머지 횟수는 Even에서 jw-Axis에 존재하는 해의 갯수
		non_JW_at_P = signChanges;
		jw_at_P = polesLength - (evenLocale - 1) - signChanges * 2;

	else
		% Even 항이 없는 경우

		% 최초의 부호값 부여
		if routhValue(1, 1) > 0
			prevSign = 1;
		else
			prevSign = -1;
		end

		% 끝까지 진행
		for i = 2:(polesLength + 1)
			if routhValue(i, 1) > 0 && prevSign == -1
				% 음수에서 양수로 변한 경우
				signChanges = signChanges + 1;
				prevSign = 1;
			elseif routhValue(i, 1) < 0 && prevSign == 1
				% 양수에서 음수로 변한 경우
				signChanges = signChanges + 1;
				prevSign = -1;
			end
		end

		% 부호의 변화 횟수는 Odd에서 RHP에 존재하는 해의 갯수
		% 위의 반복 횟수에서 나머지 횟수는 Odd에서 LHP에 존재하는 해의 갯수
		rhp_at_Q = signChanges;
		lhp_at_Q = polesLength - signChanges;
	end


	% 결과 출력
	fprintf("P_LHP: %d, P_RHP: %d, P_jw_Axis: %d\n", non_JW_at_P, non_JW_at_P, jw_at_P);
	fprintf("Q_LHP: %d, Q_RHP: %d, Q_jw_Axis: %d\n", lhp_at_Q, rhp_at_Q, 0);
	fprintf("Total_LHP: %d, Total_RHP: %d, Total_jw_Axis: %d\n", non_JW_at_P + lhp_at_Q, non_JW_at_P + rhp_at_Q, jw_at_P + 0);

	% System 특성 출력
	if (non_JW_at_P + rhp_at_Q) > 0
		fprintf("Unstable System\n");
	elseif (jw_at_P + 0) > 0
		fprintf("Marginally Stable System\n");
	else
		fprintf("Stable System\n");
	end


	% Table 반환
	outputArray = routhValue;
end
