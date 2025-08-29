function newN = Add_Crack_Node(N,sm,maxST,endL,clearL,oneCk)

% N(1,:) seg, segment no.  
% N(2,:) x, distanct
% N(3,:) k, initial stiffness  (shear flow-slip)
% N(4,:) Ef, initial stiffness  (shear flow-slip)
% N(5,:) g, < 0 for no spring and >=0 for having spring
% the values in the column (:) are the nodal values
% N(6,:) mtxST, matrix strength 

xEndClear = endL + clearL;
mtxST     = N(6,:);
xCRK      = [];
if oneCk==1
    % allow one critical crack to activate in one segment
    nSeg  = N(1,end);    % find how many segments are there
    segID = 1;
    maxOverSS = 0;
    xCrit = [];
    for i=1:length(mtxST)
        OverSress = sm(i) - mtxST(i);
        if nSeg==1
            if OverSress >= maxOverSS
                if and(N(2,i)>=endL , N(2,i)<=xEndClear)
                    maxOverSS = OverSress;
                    xCRK = N(2,i);
                end
            end
        else
            if OverSress >= maxOverSS
                if and(N(2,i)>=endL , N(2,i)<=xEndClear)
                    maxOverSS = OverSress;
                    xCrit = N(2,i);
                end
            end
            if N(1,i) > segID
                xCRK      = [xCRK,xCrit];
                segID     = segID + 1;
                maxOverSS = 0;
                xCrit     = [];
            end
            % at the last segment
            if i==length(mtxST)
                xCRK      = [xCRK,xCrit];
            end
        end
    end
else
    % allow multiple cracks to form only in the clear span, not in the grip
    for i=1:length(mtxST)
        if sm(i)>=mtxST(i)
            if and(N(2,i)>= endL , N(2,i)<= xEndClear)
                xCRK = [xCRK,N(2,i)];
            end
        end
    end
end


% add additional crack nodes in to the existing node array
newN  = N;
for i=1:length(xCRK)
    for j=1:size(N,2)
        if N(2,j) >= xCRK(i)
            newN = [N(:,1:j),N(:,j:end)];
            newN(1,j+1:end) = N(1,j:end)+1;
            newN(6,j:j+1)   = maxST;    % assign to max matrix strength for graphical purpuse in figure 4
            break;
        end
    end
    N = newN;
end
