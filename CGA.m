function optimal_parameters = CGA(FitFun,initial_guess,lower_bound,uper_bound)
%Custom Genetic Algorithm. Developed by Edwin Juárez on 2016-03-21

% tic
% disp('********************************************************************************')
% fprintf('The current date and time is\n');
% disp(fix(clock))


%Initialize the parameters
parameters = initial_guess;

chromosome_length = length(parameters);
fitness = -inf;

% fprintf('The initial guess is:');
% parameters
fitness = FitFun(initial_guess)

%  Choose parameters:
%%% { Population Size, N : Depends on the dimensions of the sample space
PopSize = 2E2;

%%% { Number of mating individuals must be an even number
FittestNumber = min(PopSize*0.1,10^3*0.05); % 10 percent of the population will reproduce up to 50 individuals

%%% { Number of "Elite" individuals who will remain from Gen[i] to Gen[i+1]
Elite = min(PopSize*0.1,10^3*0.05); %10 percent of the population up to 50 individuals

%%% { Number of Generations to Simulate: How many iterations to simulate? Presumably the more the better.
LastGen = 75;

%%% { Mutation Rates: Probability of each gene (parameter) having a point mutation.
MutRate = 0.98; % there is a 98 percent chance of a mutation. There is a lot of genetic variation!
MutMagnitude = 2;%with a mutation the parameter will change up to 100%

%%% { Crossover Points: Location(s) where the individuals swap genes in producing children (next generation).
% CrossPoint = 1; %swap (after) the first "chromosome" [parameter]
CrossPoint = floor(1+chromosome_length*rand()); % Discrete Uniform RV ranging from 1 to chromosome_length (the number of parameters)

%  Initialize the G[0] population randomly: Create a set of N solutions randomly
% Gi = RandomParameters(PopSize,parameters);
Gi = RandomParameters(PopSize,parameters,1,lower_bound,uper_bound);
vanguardia=Gi(1:4,:);
prev_vanguardia = Gi(1,:);


% figure(1)
% set(1,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); 
% hold on

%This fileID will be used store the fittest individual of each Generation
% % % % fileID_fit = fopen(strcat(pwd,'\temp_outputs\',num2str(LastGen),'_Generations_Fittness_Evolution_',line,'.txt'),'at');
%This fileID will be used store the entire Generation!
% fileID_all = fopen(strcat(pwd,'\temp_outputs\',num2str(LastGen),'_Generations_All_Data_',line,'.txt'),'at');
% diary off
i=0;
fprintf('Runing the custom genetic algorithm\n')
while(i<LastGen)

%  Calculate fitness for population
% if(i==0) fprintf('Fitness\n'); toc; end
Gi_fit = GenerationFitness(FitFun,Gi);

%  Select mates to create children for the G1 (N+1) population
%%% { Mate selection: Individuals ranked proportional to their fittness
% if(i==0) fprintf('Rank\n'); end
% Gi_fit = RankGen(Gi_fit); %Order them from the most fit to the least fit
temp = sortrows([Gi Gi_fit],-(chromosome_length+1));
Gi_fittest = temp(1:FittestNumber,1:chromosome_length); %Consider only the fittest individuals
%%% { Randomly assign mates
Gi_mate = Gi_fittest(randperm(FittestNumber),:);
Gi_mate_1 = Gi_mate(1:FittestNumber/2,:);
Gi_mate_2 = Gi_mate(FittestNumber/2+1:end,:);

%%% { Mate: Genes are exchanged prescribed by cross-over points
% if(i==0) fprintf('Mate\n');  toc; end
Offsprings = crossover(Gi_mate_1,Gi_mate_2,CrossPoint);
% Offsprings = Gi_fittest; % Turning Crossover off!!
%%% { Introduce point mutations: 
% if(i==0) fprintf('Mutate\n');  toc; end
Offsprings = mutate(Offsprings,MutRate,MutMagnitude,i);
Offsprings = CheckBoundaries(Offsprings,lower_bound,uper_bound);
%%% { Clone the Elite members and mutate the clones
Clones = mutate(Gi_fittest(1:Elite,:),MutRate,MutMagnitude,1);
Clones = CheckBoundaries(Clones,lower_bound,uper_bound);
% Clones_2 = mutate(Gi_fittest(1:round(Elite/10),:),MutRate,MutMagnitude*0.1,1);
% Clones_2 = CheckBoundaries(Clones_2,lower_bound,uper_bound);
% Clones_3 = mutate(Gi_fittest(1:round(Elite/10),:),MutRate,MutMagnitude*0.01,1);
% Clones_3 = CheckBoundaries(Clones_3,lower_bound,uper_bound);

if not(isequal(Gi_fittest(1,:),prev_vanguardia))
    prev_vanguardia = Gi_fittest(1,:);
    vanguardia = round(LinearSearch3D(FitFun,Gi_fittest(1,:),lower_bound,uper_bound));
end

%%% "Elite" fittest individuals mate with the next generation,
%%% a mutated clone of some them also carries on.

%%% on each generation a number of random individuals show up equal to the
%%% number of Elite individuals
Gi = [vanguardia; Gi_fittest(1:Elite,:); Clones; Offsprings; RandomParameters(Elite,parameters,1,lower_bound,uper_bound)];
% Gi = [vanguardia; Gi_fittest(1:Elite,:); Clones; Clones_2; Clones_3; Offsprings; RandomParameters(Elite,parameters,1,lower_bound,uper_bound)];
% Gi = [Gi_fittest(1:Elite) Clones Offsprings];
%length of Gi = FittestNumber + 3*Elite;

%plot the Current Generation's fittest individuals
% plot(i+ones(1,length(Gi_fittest)),temp(1:FittestNumber,chromosome_length+1),'p')

%print to file the Current Generation's fittest individual

% % % % fprintf(fileID_fit,'%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f',Dosage,Gi_fit(1).birth_rate,Gi_fit(1).death_rate,Gi_fit(1).clearance_rate,Gi_fit(1).fitness,i);
% % % % fprintf(fileID_fit,'\n');
% for j = 1:length(Gi)
%     fprintf(fileID_all,'%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n',Dosage,Gi(j).birth_rate,Gi(j).death_rate,Gi(j).clearance_rate,Gi(j).fitness,i);
% end

fprintf('end of generation %d fittest value = %2.2g, params = [%d,%d,%d]\n',i,Gi_fit(1),Gi_fittest(1,1),Gi_fittest(1,2),Gi_fittest(1,3));
i = i+1;
%  Repeat fitness step and selection until the maximum number of generations is reached.
end
fprintf('After %d generations, the fittest value is %2.2g, params = [%d,%d,%d]\n',i,Gi_fit(1),Gi_fittest(1,1),Gi_fittest(1,2),Gi_fittest(1,3));
optimal_parameters = Gi_fittest(1,:);
end

function Parameters = RandomParameters(PopulationSize,OrigialParameters,scale,LowerBound,UpperBound)
%PopulationSize is the number of randomized sets of parameters this function generates.
%OriginalParemeters will be a first educated guess and the rest of the
%parameters will be generated around it.
%scale is the relative value of the change, scale = 1 means  new parameters
%will be roughly the same order of magnitude; scale = 0.1 means the new
%parameters will be roughly 1/10th of the original ones.


Parameters = zeros(PopulationSize,length(OrigialParameters));

for i = 1:PopulationSize
%     Parameters(i,:) = OrigialParameters;
    for j = 1:length(OrigialParameters);
        %NOTE: round() is used because all parameters must be integers in this example
        Parameters(i,j) = round(OrigialParameters(j)*(1+scale*(2*rand-1)));
        if Parameters(i,j) < LowerBound(j)
            Parameters(i,j) = LowerBound(j);
        elseif Parameters(i,j) > UpperBound(j)
            Parameters(i,j) = LowerBound(j);
        end
    end
end

end

function Gi_fit = GenerationFitness(FitFun,Generation)
%length of Gi - Generation i;
n = length(Generation);
Gi_fit = zeros(n,1);

for i = 1:n;
    %Compute fitness of specimen i
    Gi_fit(i) = FitFun(Generation(i,:));
    if isnan(Gi_fit(i))
        Gi_fit(i) = -inf;
    end
end
end


function Offspring = crossover(Gi_mate_1,Gi_mate_2,CrossPoint)

temp = size(Gi_mate_1);
n = temp(1);
number_of_chromosomes = temp(2);

for i = 1:n
    for j = 1:CrossPoint
        Offspring(i,j)   = Gi_mate_1(i,j);
        Offspring(n+i,j) = Gi_mate_2(i,j);
    end
    for j = CrossPoint+1:number_of_chromosomes
        Offspring(i,j)   = Gi_mate_2(i,j);
        Offspring(n+i,j) = Gi_mate_1(i,j);
    end
end

end

function G_mutated = mutate(Gi,MutRate,MutMagnitude,Mutation_dampering)

G_mutated = Gi;
temp = size(Gi);
n = temp(1);
number_of_chromosomes = temp(2);
decaying_rate = 0.9;

for i = 1:n
    for j = 1:number_of_chromosomes
        if ( binornd(1,MutRate) == 1 )
%             G_mutated(i).( Gfields{j} )   = Gi(i).( Gfields{j} )*(1+MutMagnitude*(2*rand()-1));
            %%Change the mutation here!!! = Original_value*(1+Maximum_Mutation_Percentage*Unif(-1,1)*decaying_rate^Mutation_dampering) 
            %NOTE: round() is used because all parameters must be integers in this example
            G_mutated(i,j) = round((Gi(i,j) + eps )*(1+MutMagnitude*(2*rand()-1))*decaying_rate^Mutation_dampering);
        else
            G_mutated(i,j) = Gi(i,j);
        end
    end
end

end

function Gi = CheckBoundaries(Gi,lower_bound,uper_bound)
temp = size(Gi);
n = temp(1);
number_of_chromosomes = temp(2);

for i = 1:n
    for j = 1:number_of_chromosomes
        if ( Gi(i,j)<lower_bound(j) )
            Gi(i,j) = lower_bound(j);
        end
        if( Gi(i,j)>uper_bound(j) )
            Gi(i,j) = uper_bound(j);
        end
    end
end

end

function vanguardia = LinearSearch3D(FitFun,Gi,lb,ub)

temp = size(Gi);
n = temp(1);
number_of_chromosomes = temp(2);

percent = 0.1;

vanguardia = zeros(n*4,number_of_chromosomes);

for i = 1:n
    % search the best x:
%     disp('searching x')
    vanguardia((i-1)*4+1,:) = LinearSearch(FitFun,Gi(i,:),1,lb*(1-percent),ub*(1+percent));
    % search the best y:
%     disp('searching y')
    vanguardia((i-1)*4+2,:) = LinearSearch(FitFun,Gi(i,:),2,lb*(1-percent),ub*(1+percent));
    % search the best z:
%     disp('searching z')
    vanguardia((i-1)*4+3,:) = LinearSearch(FitFun,Gi(i,:),3,lb*(1-percent),ub*(1+percent));
    % Do a local search on nearby
%     disp('searching local ball')
    vanguardia((i-1)*4+4,:) = LocalSearch(FitFun,Gi(i,:),percent*0.1,lb,ub);
%     disp('done searching local ball')
end

end

function max_x = LocalSearch(FitFun,initial_guess,percent,lb,ub)

lba = CheckBoundaries(round(initial_guess*(1-percent)),lb,ub);
uba = CheckBoundaries(round(initial_guess*(1+percent)),lb,ub);

counter = 1;
max_fit = -inf;
max_x = [0,0,0];
for a = lba(1):uba(1)
    for b = lba(2):uba(2)
        for c = lba(3):uba(3)
            temp_x = [a,b,c];
            temp_fit = FitFun(temp_x);
            if(temp_fit>max_fit)
                max_fit = temp_fit;
                max_x = temp_x;
            end
            counter = counter + 1;
        end
    end
end

end

function max_x = LinearSearch(FitFun,initial_guess,index,lb,ub)

% counter = 1;
max_fit = -inf;
max_x = initial_guess;
temp_x = initial_guess;
for variable_value = lb(index):ub(index)
    temp_x(index) = variable_value;
    temp_fit = FitFun(temp_x);
    if(temp_fit>max_fit)
        max_fit = temp_fit;
        max_x = temp_x;
    end
%     counter = counter + 1;
end

end