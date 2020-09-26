function [A] = MxMake_1050044(mx_id,n)
str1="had";
str2="trihad";
str3="toep";
str4="mc";
str5="wathen";
str6="CollegeMsg";

%mitrwo Hadamard
if strcmp(mx_id,str1)
    A=hadamard(n);
    disp(A);
    
    return
end
%trihad

if strcmp(mx_id,str2)
    A=hadamard(n);
    A=triu(A);
    disp(A);
    return
end
%toep

if strcmp(mx_id,str3)
    A=toeplitz([4,-1,zeros(1,n-2)]);
    disp(A);
    return
end
%mc
if strcmp(mx_id,str4)
    thita=1;
    s=2;
    for i=1:n
        fill=1+(i^thita);
        A(n,n)=fill;
    end
    for i=1:n
        for j=1:n
            if(i~=j)
                fill=1/(abs(i-j)^s);
                A(i,j)=fill;
            end
            if(i==j)
                fill=1+(i^thita);
                A(i,j)=fill;
            end
        end
    end
    disp(A);
    return
end
%wathen
if strcmp(mx_id,str5)
    m=11;
    A=gallery('wathen',n,m);
    disp(A);
    return
end

%collegeMsg
if strcmp(mx_id,str6)
    struct=load('CollegeMsg.mat');
    %C=(struct2cell(A.Problem.A));
    A=struct.Problem.A;
    return
end
end














