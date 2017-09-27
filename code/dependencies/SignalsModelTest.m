% Test suite for SignalsModel

% Setup
signals = SignalsValues;
signals.addSignal('a', 1:10);
signals.addSignal('b', 101:110);
signals.addSignal('c', 1001:1010);
signals.truncate(10);

mdl = SignalsModel(signals);
mdl.setInputs({'a', {'b', [2 1 0]}, {'c', [2 1]}});
mdl.setOutput('c');

%% Test simple IO vectors
[X, Y] = mdl.getIOVectors();
X0 = [(1:10)', [NaN, NaN, 101:108]', [NaN, 101:109]', (101:110)', [NaN, NaN, 1001:1008]', [NaN, 1001:1009]'];
Y0 = (1001:1010)';
assert(isequaln(X, X0) && isequaln(Y, Y0));

%% Test IO vectors with NaNs removed
[X, Y] = mdl.getIOVectors([], 'removeNaNs', true);
X0 = [(3:10)', (101:108)', (102:109)', (103:110)', (1001:1008)', (1002:1009)'];
Y0 = (1003:1010)';
assert(isequaln(X, X0) && isequaln(Y, Y0));

%% Test IO vectors with indices and stepsahead
[X, Y] = mdl.getIOVectors(5, 'stepsahead', 2);
X0 = [7, 105, 106, 107, 1005, 1006];
Y0 = 1007;
assert(isequaln(X, X0) && isequaln(Y, Y0));

[X, Y] = mdl.getIOVectors([5 6], 'stepsahead', 2);
X0 = [7, 105, 106, 107, 1005, 1006;
    8, 106, 107, 108, 1006, 1007];
Y0 = [1007; 1008];
assert(isequaln(X, X0) && isequaln(Y, Y0));

%% Test IO vectors with substitutes
inputsubs = struct;
inputsubs.b.start = 6;
inputsubs.b.values = [206; 207];
inputsubs.c.start = 5;
inputsubs.c.values = 2005;
[X, Y] = mdl.getIOVectors(5, 'stepsahead', 2, 'inputsubs', inputsubs);
X0 = [7, 105, 206, 207, 2005, 1006];
Y0 = 1007;
assert(isequaln(X, X0) && isequaln(Y, Y0));

%% Test IO vectors with substitutes and Symbolic toolbox
if license('test', 'Symbolic_Toolbox')
    inputsubs = struct;
    inputsubs.b.start = 6;
    inputsubs.b.values = sym('b', [2, 1]);
    inputsubs.c.start = 5;
    inputsubs.c.values = sym('c', [3,1]);
    [X, Y] = mdl.getIOVectors(5, 'stepsahead', 2, 'inputsubs', inputsubs, 'type', 'sym');
    X0 = [7, 105, inputsubs.b.values(1), inputsubs.b.values(2), inputsubs.c.values(1), inputsubs.c.values(2)];
    Y0 = inputsubs.c.values(3);
    assert(isequaln(X, X0) && isequaln(Y, Y0));
end

%% Test IO vectors with substitutes and CasADi
if exist('casadi.SX', 'class')
    inputsubs = struct;
    inputsubs.b.start = 6;
    inputsubs.b.values = casadi.SX.sym('b', [2, 1]);
    inputsubs.c.start = 5;
    inputsubs.c.values = casadi.SX.sym('c', [3,1]);
    [X, Y] = mdl.getIOVectors(5, 'stepsahead', 2, 'inputsubs', inputsubs, 'type', 'casadi.SX');
    X0 = [7, 105, inputsubs.b.values(1), inputsubs.b.values(2), inputsubs.c.values(1), inputsubs.c.values(2)];
    Y0 = inputsubs.c.values(3);
    assert(X.is_equal(X0) && Y.is_equal(Y0));
end