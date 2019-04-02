function lab = HGxyz2lab(xyz,wp)

    X = xyz(:,1)/wp(1);
    Y = xyz(:,2)/wp(2);
    Z = xyz(:,3)/wp(3);

    fX = real(X.^(1/3));
    i = (X < 0.008856);
    fX(i) = X(i)*(841/108) + (4/29);

    fY = real(Y.^(1/3));
    i = (Y < 0.008856);
    fY(i) = Y(i)*(841/108) + (4/29);

    fZ = real(Z.^(1/3));
    i = (Z < 0.008856);
    fZ(i) = Z(i)*(841/108) + (4/29);

    lab = zeros(size(xyz));
    lab(:,1) = 116*fY - 16;    
    lab(:,2) = 500*(fX - fY);  
    lab(:,3) = 200*(fY - fZ);  

end
