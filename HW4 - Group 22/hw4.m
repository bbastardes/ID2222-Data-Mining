E = csvread('data/example1.dat');
k = 4;

%E = csvread('example2.dat');
%k = 2;

sigma = 1;

% Adjacency matrix
col1 = E(:,1);
col2 = E(:,2);
max_ids = max(max(col1,col2));
%As= sparse(col1, col2, 1, max_ids, max_ids); 
%A = full(As)

G = graph(col1, col2);
A = adjacency(G);
A = full(A);


%for i = 1:241
%    for j = 1:241
%        A(i,j) = exp(-norm((i-j),2)/2*sigma^2);
%        
%    end
%end

% Diagonal matrix and L
D = diag(sum(A,2));
L = D^(-1/2)*A*D^(-1/2);

% k largest eigenvectors of L and matrix X
[evec, eval] = eig(L);
[eval_sorted, idx] = sort(diag(eval),'descend');
X = evec(:,idx(1:k));

% Y matrix normalising X
Y = bsxfun(@rdivide, X, sqrt(sum(X.^2)));
Y(isnan(Y)) = 0;

% clustering
color = [1 0 0; 0 1 0; 0 0 1; 1 0 1];
[idx,C] = kmeans(Y,k);
G = graph(col1,col2);
figure,
h = plot(G,'NodeColor',color(idx,:));


