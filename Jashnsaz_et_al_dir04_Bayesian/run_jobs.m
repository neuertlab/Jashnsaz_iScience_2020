function run_jobs()
% steps 10 chains job_ID: 1-10
% diverse 10 chains job_ID: 11-20
% D-Optimal 10 chains job_ID: 21-30
job_ID = 1; 
%for job_ID=1:30 % all 30 jobs
    run_Bayes(job_ID)
%end
end