function [x] = SMW_solve_1050044(A,b,M,P,Q,sdir)
str1="rowwise";
str2="colwise";
% if(sdir==str1)||(sdir==str2)%tha ekteleitai mono gia rowwise kai colwise
%     %A = randn(10,10);           %gia opoiadipote alli timi sdir tha pairnei tis eisodous tis SMW_solve
%     %b = randn(10,1);

%    %C==M
%     A0=sparse(A0);%Αραιό μητρώο μόνο με μη μηδενικά στοιχεία
%     %για εξοικονόμηση πράξεων
%
x = M \ b;%lisi sistimatos Mx=b gia to x kai vima 1
y = M \ P;%lisi sistimatos My=A0 gia to y kai vima 2
%     v = eye(n)%Q
%     v=sparse(v);%Αραιό μητρώο μόνο με μη μηδενικά στοιχεία
%     %για εξοικονόμηση πράξεων
% end
m = size(A);
n = m(1);
%rowwise
v=Q;
if sdir==str1
    U=eye(n);%katalliles energeies gia to  P kai Q
    v=U';
    %v=sparse(v);
    for l=1:n-1
        x=x-(v(:,l)'*x/(1+v(:,l)'*y(:,l)))*y(:,l);
        for k=(l+1):n
            y(:,k)=y(:,k)-(v(:,l)'*y(:,k)/(1+v(:,l)'*y(:,l)))*y(:,l);
        end
    end
    
    x=x-(v(:,n)'*x/(1+v(:,n)'* y(:,n)))*y(:,n);
    return
end

%colwise
if sdir==str2
    for l=1:n-1
        x=x-(v(:,l)'*x/(1+v(:,l)'*y(:,l)))*y(:,l);
        for k=(l+1):n
            y(:,k)=y(:,k)-(v(:,l)'*y(:,k)/(1+v(:,l)'*y(:,l)))*y(:,l);
        end
    end
    
    
    x=x-(v(:,n)'*x/(1+v(:,n)'* y(:,n)))*y(:,n);
    return
end
%auti einai apantisi dia to "diaforetika" otan den einai colwise I rowwise
if(sdir~=str1)||(sdir~=str2)
    %A0 = diag(diag(A));
    m = size(Q);
    n = m(1);
    v = eye(n);
    %v=sparse(v);
    %M =sparse(M);
    %A = M + Q*v';%Αντικατάσταση του Α με επίλυση Μ+PQ'
    %x  =  A \ b;
    for l=1:n-1
        x=x-(v(:,l)'*x/(1+v(:,l)'*y(:,l)))*y(:,l);
        for k=(l+1):n
            y(:,k)=y(:,k)-(v(:,l)'*y(:,k)/(1+v(:,l)'*y(:,l)))*y(:,l);
        end
    end
    
    
    x=x-(v(:,n)'*x/(1+v(:,n)'* y(:,n)))*y(:,n);
    return
    
end

%i Η λογική μου είναι ότι ο υπολογισμός του inv(Μ) ειναι λιγότερο
%δαπανηρός απο ότι αν έκανα inv(A) γιατί το Μ ειναι δισδιαγώνιο του Α
%άρα ειναι πιο αραιό μητρώο απο το Α άρα λιγότεροι υπολογισμοί
%Επίσης ίσως να μπορούσα να χρησιμοποιήσω Lu παραγοντοποίηση για να
%διασπάσω τα μητρώα
%Κατόπιν πειραμάτων με την timeit και toc διαπίστωσα ότι ο υπολογισμός
%των αναστρόφων μια φορά ειναι πιο οικονομικός απο τους συνεχείς
%υπολογισμούς.Επίσης διαπίστωσα ότι ειναι πιο οικονομικό ακόμα και απο το
% να ανάτρεξει στη μνήμη για ανάκτηση του αποτελέσματος.

%COLWISE PIO PANW
%ROWISE TO U=eye(n) kai v=U' THA
%DIAFORETIKA THA KANW M+PQ'*X(AO+PQ
%THA TO VAZW KAI STO COLWISE KAI STO ROWWISE
%MONO STO ELSE DEN THA TO VAZW

end


