function t = myTrace(X, dims)
%MYTRACE  견고한 Trace 계산 함수
%
%   t = MYTRACE(X)        : 2-D 행렬 X의 주대각합(Trace)을 반환합니다.
%   t = MYTRACE(X, DIMS)  : N-D 배열 X에 대해, 행 방향 차원과 열 방향
%                           차원을 각각 DIMS(1), DIMS(2)로 지정하여 Trace
%                           를 계산합니다. (예: DIMS = [2 4])
%
% ──────────────────────────────────────────────────────────────────────────
% 특징
%   • 직사각형 행렬 허용  →  min(size(X)) 길이만큼의 주대각선 합
%   • N-D 페이지 지원     →  나머지 차원은 독립 페이지로 간주
%   • 실수·복소·논리·희소·심볼릭 모두 지원
%   • 빈 행렬 입력 시 0 반환 (MATLAB 내장 trace와 동일)
%   • 내장 trace에 이름이 가려지거나 경로 충돌이 있을 때 대체 사용
%
% 작성: 2025-07-14  (박정훈 님 요청)
% ──────────────────────────────────────────────────────────────────────────

    narginchk(1,2);                                     % 인수 확인
    validateattributes(X, {'numeric','logical','sym'}, ...
                          {'nonempty'}, mfilename, 'X', 1);

    % 기본 행·열 차원 결정
    if nargin < 2
        dims = [1 2];
    end
    assert(numel(dims)==2, 'DIMS 인수는 행·열 차원 두 개여야 합니다.');
    rowDim = dims(1);
    colDim = dims(2);
    sz     = size(X);

    % 입력 차원 범위 점검
    maxDim = max([rowDim colDim ndims(X)]);
    if numel(sz) < maxDim, sz(maxDim) = 1; end
    assert(rowDim ~= colDim, '행 차원과 열 차원이 같을 수 없습니다.');

    % 행/열 길이 및 대각선 길이
    m   = sz(rowDim);
    n   = sz(colDim);
    k   = min(m,n);                     % 실제 대각선 원소 수

    % 행·열 차원을 앞으로 옮겨서 재배열   (permute는 복사 없는 뷰)
    order = 1:maxDim;
    order([1 rowDim]) = order([rowDim 1]);   % rowDim → 1
    order([2 colDim]) = order([colDim 2]);   % colDim → 2
    Xp    = permute(X, order);               % 크기: m × n × P

    % 나머지 차원을 납작하게 묶어 페이지화
    P   = prod(sz(setdiff(1:maxDim, [rowDim colDim])));
    Xp  = reshape(Xp, m, n, P);              % 3-D: m × n × P

    % 대각선 인덱스 (희소 행렬도 안전)
    diagIdx = 1:(m+1):((k-1)*(m+1)+1);

    % 벡터화된 대각합
    if issparse(X)
        % 희소행렬은 full로 변환하지 않고 spdiags 이용
        t = sum(spdiags(Xp(:,1:k,:), 0), 1);
    else
        % 일반/심볼릭/복소
        t = sum(Xp(diagIdx, :, :), 1);       % m × n × P → 1 × n × P
    end

    % 페이지별 스칼라 결과로 변환, 원래 페이지 차원 형태 복원
    t = squeeze(t);                         % → 1 × P
    outSz = sz;
    outSz([rowDim colDim]) = [];            % 행·열 차원 제거
    t = reshape(t, outSz);
end
