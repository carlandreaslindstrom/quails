% Saves the solution provided to file
function [] = saveSolution(best_x, fileID)

    % open the solutions file
    fileID = fopen(['solutions/solutions_' fileID '.dat'],'w');
    
    % print the solution
    fprintf(fileID,'%1.18e\n', best_x);
    
    % close the file
    fclose(fileID);

end
