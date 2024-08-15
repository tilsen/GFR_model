function [D] = prep_data(D)


if isfield(D,'t') && ~isfield(D,'dt')
    D.dt = D.t{1}(2)-D.t{1}(1);
elseif ~isfield(D,'t') && isfield(D,'dt')
    D.t = {0:D.dt:((length(D.y{1})-1)*D.dt)};
end

if ~isfield(D,'pw_t0')
    D.pw_t0{1} = D.t(1);
    D.pw_t1{1} = D.t(end);
end

D.Npw = length(D.pw_t0{1});

for i=1:D.Npw
    D.y_pw{1}{i} = D.y{1}(D.t{1}>=D.pw_t0{1}(i) & D.t{1}<D.pw_t1{1}(i));
end

D = struct2table(D);


D.pw_dur{1} = D.pw_t1{1}-D.pw_t0{1};

end

