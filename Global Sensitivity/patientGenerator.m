function[fakePatients] = patientGenerator(n, lowMult, highMult, fcr1, fcr2, disease, bodycell)

[p, ~, ~, ~] = MultipleFcRParameters(fcr1, fcr2, disease, bodycell);

% p(9) = mean([35e-6 52e-6]);
% p(11) = mean([1e-6 4.5e-6]);
% p(13) = mean([9.1e-6 8.9e-6]);
% p(15) = mean([2.1e-6 1.7e-6]);

pRange = [p * lowMult; p * highMult];

pOptions = zeros(n,31);

for i = 1:31
    pmin = pRange(1,i);
    pmax = pRange(2,i);
    pOptions(:,i) = logspace(log10(pmin),log10(pmax),n);
end

fakePatients = zeros(n,31);

for i = 1:31
    fakePatients(:,i) = randperm(n)';
end

for i = 1:n
    for j = 1:31
        fakePatients(i,j) = pOptions(fakePatients(i,j),j);
    end
end