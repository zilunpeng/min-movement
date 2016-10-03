prompt = 'enter number of pts ';
num_pts = input(prompt);
pts_coord = zeros(num_pts,2);

for i=1:num_pts
    prompt = 'enter x coordinate ';
    pts_coord(i,1) = input(prompt);
end

prompt = 'enter delta ';
delta = input(prompt);

p = (1:num_pts);
all_possible_order = perms(p);
% all_possible_order = [1 3 5 2 4 6];
[num_possible_way, num_pts] = size(all_possible_order);

lifted_pts_coord = zeros(num_pts,2);
best_lifted_pts_coord = zeros(num_pts,2);
current_total_movement = 0;
min_total_movement = -log(0);

%find best soln
for i=1:num_possible_way
    lifted_pts_coord(1,1) = pts_coord(all_possible_order(i,1),1);
    for j=2:num_pts
        %get x-coord of current point
        x = pts_coord(all_possible_order(i,j),1);
        %find conlifts with pts that have been considered
        [conflict_pts_r_ind, conflict_pts_c_ind] = find(or(lifted_pts_coord(1:j-1,1) > x & lifted_pts_coord(1:j-1,1) < x+delta, lifted_pts_coord(1:j-1,1) < x & lifted_pts_coord(1:j-1,1) > x-delta));
        %if current pt has no conflict, simply add it; otherwise, find
        %conflict pt with max y coord and put current pt above it
        if isempty(conflict_pts_r_ind) == 1
            lifted_pts_coord(j,1) = x;
        else
            %multiply by 2 to get column indicies of conflict pts' y coords
            conflict_pts_c_ind = conflict_pts_c_ind * 2;
            %access y_coords of conflict pts
            ind = sub2ind(size(lifted_pts_coord), conflict_pts_r_ind, conflict_pts_c_ind);
            conflict_pt_y_coord = lifted_pts_coord(ind);
            %access x_coord of conflict pts
            ind = sub2ind(size(lifted_pts_coord), conflict_pts_r_ind, conflict_pts_c_ind/2);
            conflict_pt_x_coord = lifted_pts_coord(ind);
            %calculate the lifted distances to avoid all conflict pts
            num_conflict_pts = length(conflict_pt_x_coord);
            y_lifted = -log(0);
            for k=1:num_conflict_pts
               y_temp_lifted = sqrt(delta^2 - (x-conflict_pt_x_coord(k))^2)+conflict_pt_y_coord(k);
               dists_to_conflict_pts = sqrt((ones(num_conflict_pts,1)*x-conflict_pt_x_coord).^2 + (ones(num_conflict_pts,1)*y_temp_lifted-conflict_pt_y_coord).^2);
               if isempty(find(round(dists_to_conflict_pts,3) < delta, 1)) && y_temp_lifted<y_lifted
                   y_lifted = y_temp_lifted;
               end                
            end
%             [max_conflict_y_lifted, ind] = max(conflict_pt_y_coord);
%             if length(find(conflict_pt_y_coord == max_conflict_y_lifted)) > 1
%                 dists_to_x = abs(conflict_pt_x_coord(conflict_pt_y_coord == max_conflict_y_lifted) - x);
%                 y_lifted = sqrt(delta^2 - (min(dists_to_x)^2)) + max_conflict_y_lifted;
%             else
%                 y_lifted = sqrt(delta^2 - (conflict_pt_x_coord(ind)-x)^2) + max_conflict_y_lifted;
%             end
            lifted_pts_coord(j,1) = x;
            lifted_pts_coord(j,2) = y_lifted;
            current_total_movement = current_total_movement + y_lifted;
        end
    end
    if current_total_movement < min_total_movement
        min_total_movement = current_total_movement;
        best_lifted_pts_coord = lifted_pts_coord;
    end
    current_total_movement = 0;
    lifted_pts_coord = zeros(num_pts,2);
end


%plot
for i=1:num_pts
    drawCircle(best_lifted_pts_coord(i,1),best_lifted_pts_coord(i,2),delta/2);
    line([best_lifted_pts_coord(i,1) best_lifted_pts_coord(i,1)],[-max(best_lifted_pts_coord(:,2)) max(best_lifted_pts_coord(:,2))+2*delta])
    text(best_lifted_pts_coord(i,1),best_lifted_pts_coord(i,2),num2str(best_lifted_pts_coord(i,2)))
    text(best_lifted_pts_coord(i,1),-delta,num2str(best_lifted_pts_coord(i,1)))
    hold on
end
axis([0-delta max(best_lifted_pts_coord(:,1))+delta 0-delta max(best_lifted_pts_coord(:,2))+delta])
