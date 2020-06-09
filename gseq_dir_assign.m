function output=gseq_dir_assign(arraysize, flag)
    n=(arraysize+1)/2;
    arraysize=2*n-1;
    sequence=zeros(2,arraysize^2);
    sequence(1,1)=n;
    sequence(2,1)=n;
    dx=+1;
    dy=-1;
    stepx=+1;
    stepy=-1;

    direction=+1; %+1 means x direction, -1 means y direction
    counter=0;
    for i=2:arraysize^2
        counter=counter+1;

        if direction==+1
            sequence(1,i)=sequence(1,i-1)+dx;
            sequence(2,i)=sequence(2,i-1);
            if counter==abs(stepx)
                counter=0;
                direction=direction*-1;
                dx=dx*-1;
                stepx=stepx*-1;
                if stepx>0
                    stepx=stepx+1;
                else
                    stepx=stepx-1;
                end

            end

        else
            sequence(1,i)=sequence(1,i-1);
            sequence(2,i)=sequence(2,i-1)+dy;

            if counter==abs(stepy)
                counter=0;
                direction=direction*-1;
                dy=dy*-1;
                stepy=stepy*-1;

                if stepy>0
                    stepy=stepy+1;
                else
                    stepy=stepy-1;
                end

            end

        end
    end
    seq=(sequence(1,:)-1)*arraysize+sequence(2,:);
    seqf(1,1:arraysize^2)=seq;   
    output = seqf;
    as = (1:2:arraysize).^2;
    MOVE_NUM = zeros([1, length(as)-1]);
    if flag == 'N'
        MOVE_NUM = (as(2:end)-as(1:(end-1)))/2;
    elseif flag == 'NE'
        MOVE_NUM = 3*(as(2:end)-as(1:(end-1)))/8;
    elseif flag == 'E'
        MOVE_NUM = (as(2:end)-as(1:(end-1)))/4;
    elseif flag == 'SE'
        MOVE_NUM = (as(2:end)-as(1:(end-1)))/8;
    elseif flag == 'S'
        
    elseif flag == 'SW'
        MOVE_NUM = 7*(as(2:end)-as(1:(end-1)))/8;
    elseif flag == 'W'
        MOVE_NUM = 3*(as(2:end)-as(1:(end-1)))/4;
    elseif flag == 'NW'
        MOVE_NUM = 5*(as(2:end)-as(1:(end-1)))/8;
    end
    for i = 2:length(as)
        output((as(i-1)+1):as(i)) = circshift(seqf((as(i-1)+1):as(i)), [0, MOVE_NUM(i-1)]);
    end
end

