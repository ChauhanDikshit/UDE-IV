function X = insert(A, b, pos)

X = [A(1:pos-1,:); b; A(pos:end-1,:)];

end

