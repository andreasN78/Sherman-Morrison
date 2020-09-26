clc;
clear all;
prompt = {'Δώστε το κατάλληλο μέγεθος πίνακα:','Δώστε τον τύπο του μητρώου:','Δώστε το sdir:'};
title = 'Είσοδος Δεδομένων';
dims = [1 35];
definput = {'4','had','colwise'};
answer = inputdlg(prompt,title,dims,definput);
%answer=cell2table(answer)
n=answer(1);
sdir=answer{3};
n=str2double(n);
mx_id=answer(2);
A=MxMake_1050044(mx_id,n);

str1='rowwise';
str2='colwise';
str3='CollegeMsg';

%A*x=b
nsize=size(A);
n=nsize(2);


xsol=ones(1,n);
for k=1:n/2
    thesis1=(2*k)-1;
    xsol(thesis1)=1;
end
for k=1:n/2
    thesis2=2*k;
    xsol(thesis2)=(((-1)^(k+1))*(1/(2*k)));
end

xsol=xsol';
disp(A);
disp(xsol);
b = A*xsol;
x_mat = A \ b;%Το αποτέλεσμα που πρέπει να βρώ
%A = randn(10,10);           %gia opoiadipote alli timi sdir tha pairnei tis eisodous tis SMW_solve
m = size(A);
nA = m(1);
M= diag(diag(A));%C==M
% M=sparse(M);%Αραιό μητρώο μόνο με μη μηδενικά στοιχεία
%για εξοικονόμηση πράξεων
P = A -M ;%P
x = M \ b;%lisi sistimatos Mx=b gia to x kai vima 1
y = M \ P;%lisi sistimatos My=A0 gia to y kai vima 2
Q = eye(nA);%Q
%Q=sparse(Q);%Αραιό μητρώο μόνο με μη μηδενικά στοιχεία
%για εξοικονόμηση πράξεων
switch sdir
    case {'colwise'}
        %nothing
    case {'rowwise'}
        %nothing
    otherwise
        m = size(Q);
        nA = m(1);
        v = eye(nA);
        %v=sparse(v);
        
        A = M + Q*v';%Αντικατάσταση του Α με επίλυση Μ+PQ'
end

xSmw=SMW_solve_1050044(A,b,M,P,Q,sdir)%Το αποτέλεσμα του Αλγόριθμου 1
%Εξακρίβωση αποτελέσματος
akrivwsEmprosSfalma=norm(x-xsol)/norm(xsol)
%sfalma=(xSmw-xsol)/norm(xSmw);
%δείκτης κατάστασης για μητρώο Α
deiktisK=condest(A)
%α posteriori Πίσω σφάλμα
aposteriori=norm(b-(A*x_mat))/(norm(A,'fro')*norm(x_mat,'fro'))+norm(b)
%εμπρός φράγμα
emprosfragma=(2*aposteriori*deiktisK)

%Εξακρίβωση αποτελέσματος
akrivwsEmprosSfalmaSmw=norm(x-xsol)/norm(xsol)
%sfalma=(xSmw-xsol)/norm(xSmw);

%α posteriori Πίσω σφάλμα
aposterioriSmw=norm(b-(A*xSmw),'fro')/(norm(A,'fro')*norm(xSmw,'fro'))+norm(b)
%εμπρός φράγμα
emprosfragmaSmw=(2*aposteriori*deiktisK)
