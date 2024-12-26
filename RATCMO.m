classdef RATCMO < ALGORITHM
% <multi> <real> <constrained>

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population{1} = Problem.Initialization(); % Main task
            Fitness{1}    = CalFitness(Population{1}.objs,Population{1}.cons);
            Population{2} = Problem.Initialization(); % Global auxiliary task
            Fitness{2}   = CalFitness(Population{2}.objs);
            Population{3} = Problem.Initialization(); % Local auxiliary task
            Fitness{3}   = CalFitness(Population{3}.objs,Population{3}.cons);
            %% Evaluate the Population
            RMP = 0.8;
            Nt = 2;      % The number of tasks
            Beta = 0.3;   % Learning phase
            Pmin = 0.1;   % Minimum selection probability
            A = Population{1};
            rwd   =  zeros(Nt, 1);                 % Rewards

            cons = Population{1}.cons;
            cons(cons<0) = 0;
            VAR0 = max(sum(cons,2));
            if VAR0 == 0
                VAR0 = 1;
            end
            X=0;
            %% Optimization
            while Algorithm.NotTerminated(A)
                Q = [];
                cp  = (-log(VAR0)-6)/log(1-0.5);
                VAR = VAR0*(1-X)^cp;
                % Calculate the selection probability
                if Problem.FE <= Beta * Problem.maxFE
                    % Stage 1: Evolution stage
                    pro  =  1 / Nt * ones(1, Nt);
                else
                    % Stage 2: Competition stage
                    if sum(rwd) ~= 0
                        pro   =   Pmin / Nt + (1 - Pmin) * rwd ./ sum(rwd);
                        pro   =   pro ./ sum(pro);
                    else
                        pro   =   1 / Nt * ones(1, Nt);
                    end
                end
                % Determine the a task based on the selection probability using roulette wheel method
                r = rand;
                for t = 1:Nt
                    if r <= sum(pro(1:t))
                        k = t;
                        break;
                    end
                end
                %for i = 1:Problem.N
                % Choose the parents
                %P = randperm(Problem.N);
                MatingPool1 = TournamentSelection(2,Problem.N/2,Fitness{1});
                %MatingPool11 = TournamentSelection(2,Problem.N/2,Fitness{1});
                Offspring1(1:Problem.N/4) = OperatorGAhalf(Problem,[Population{1}(MatingPool1)]);
                Offspring1(Problem.N/4+1:Problem.N/2) = OperatorDE(Problem,Population{1}(1:4:Problem.N),Population{1}(MatingPool1(1:end/2)),Population{1}(MatingPool1(1+end/2:end)));

                if k == 1
                    if rand < RMP
                        MatingPool2 = TournamentSelection(2,Problem.N,Fitness{2});
                        Offspring2 = OperatorGAhalf(Problem,[Population{2}(MatingPool2)]);
                        Offspring=[Offspring1,Offspring2];
                    else
                        MatingPool3 = TournamentSelection(2,min(length(Population{3}),Problem.N/2),Fitness{3});
                        Offspring3 = OperatorGA(Problem,[Population{3}(MatingPool3)]);
                        Offspring=[Offspring1,Offspring3];
                    end
                    for j=1:3
                        if j == 1
                            [Population{j},Fitness{j}] = EnvironmentalSelection( [Population{1},Offspring],Problem.N,true);
                        elseif j==2
                            [Population{j},Fitness{j}] = EnvironmentalSelection( [Population{2},Offspring],Problem.N,false);
                        else
                            [Population{j},Fitness{j}] = EnvironmentalSelection2([Population{3},Offspring],Problem.N,VAR);
                        end
                    end
                end
                if k==2
                    if rand < RMP
                        MatingPool3 = TournamentSelection(2,min(length(Population{3}),Problem.N/2),Fitness{3});
                        Offspring3 = OperatorGA(Problem,[Population{3}(MatingPool3)]);
                        Offspring=[Offspring1,Offspring3];
                    else
                        MatingPool2 = TournamentSelection(2,Problem.N,Fitness{2});
                        Offspring2 = OperatorGAhalf(Problem,[Population{2}(MatingPool2)]);
                        Offspring=[Offspring1,Offspring2];
                    end
                    for j=1:3
                        if j == 1
                            [Population{j},Fitness{j}] = EnvironmentalSelection( [Population{1},Offspring],Problem.N,true);
                        elseif j==2
                            [Population{j},Fitness{j}] = EnvironmentalSelection( [Population{2},Offspring],Problem.N,false);
                        else
                            [Population{j},Fitness{j}] = EnvironmentalSelection2([Population{3},Offspring],Problem.N,VAR);
                        end
                    end
                end
                num = size(Offspring1,2);

                Q = [Q Offspring];
                %end
                s = size(A,2);
                [A,Next] =  ArchiveUpdate([A Q],Problem.N);

                if size(Q,2) > 0 && Problem.FE > 0.1 * Problem.maxFE
                    % update Archive
                    if s >= Problem.N
                        rwd(k)   =  rwd(k) + sum(Next(Problem.N+num+1:end))/Problem.N;   % update the reward of the main task
                    end
                end
                X=Problem.FE/Problem.maxFE;
            end
        end
    end
end
