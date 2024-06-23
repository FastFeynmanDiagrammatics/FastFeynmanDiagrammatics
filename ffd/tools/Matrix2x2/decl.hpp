namespace ffd::matrix_2x2{
	
    class Matrix2x2 {
    public: std::array<Real, 4> coeff;
            Matrix2x2 () {};
            Matrix2x2 (Real x){coeff[0]=x;coeff[1]=0.;coeff[2]=0.;coeff[3]=x;};
            Matrix2x2 operator+(const Matrix2x2 &) const;
            Matrix2x2 operator-(const Matrix2x2 &) const;
            Matrix2x2 operator+=(const Matrix2x2 &);
            Matrix2x2 operator-=(const Matrix2x2 &);
            Matrix2x2 operator*(const Matrix2x2 &) const;
            Matrix2x2 operator*=(const Matrix2x2 &);
            Matrix2x2 operator/=(Real x);
            Real det(){return coeff[0] * coeff[3] - coeff[1] * coeff[2];};
};

Matrix2x2 Matrix2x2::operator+(const Matrix2x2 & M) const{
    Matrix2x2 temp = *this;
    for (uint i=0; i<4; ++i){
        temp.coeff[i] += M.coeff[i];  
    }
    return temp; 
}

Matrix2x2 Matrix2x2::operator-(const Matrix2x2 & M) const{
    Matrix2x2 temp = *this;
    for (uint i=0; i<4; ++i){
        temp.coeff[i] -= M.coeff[i];  
    }
    return temp; 
}

Matrix2x2 Matrix2x2::operator*(const Matrix2x2 & M) const{
    Matrix2x2 temp;
    temp.coeff[0] = this->coeff[0]*M.coeff[0] + this->coeff[1]*M.coeff[2];  
    temp.coeff[1] = this->coeff[0]*M.coeff[1] + this->coeff[1]*M.coeff[3];  
    temp.coeff[2] = this->coeff[2]*M.coeff[0] + this->coeff[3]*M.coeff[2];  
    temp.coeff[3] = this->coeff[2]*M.coeff[1] + this->coeff[3]*M.coeff[3];  
    return temp; 
}

Matrix2x2 Matrix2x2::operator+=(const Matrix2x2 &  M){
    return *this = *this + M; 
}

Matrix2x2 Matrix2x2::operator-=(const Matrix2x2 & M){
    return *this = *this - M; 
}

Matrix2x2 Matrix2x2::operator*=(const Matrix2x2 & M){
    return *this = *this * M; 
}

Matrix2x2 operator/(Real x, const Matrix2x2 & M){
   Real det = M.coeff[0] * M.coeff[3] - M.coeff[1] * M.coeff[2];
   Real prefactor = x / det;
   Matrix2x2 temp;
   temp.coeff[0] = prefactor *  M.coeff[3];
   temp.coeff[1] = prefactor * -M.coeff[1];
   temp.coeff[2] = prefactor * -M.coeff[2];
   temp.coeff[3] = prefactor *  M.coeff[0];
   return temp;
}

Matrix2x2 Matrix2x2::operator/=(Real x) {
   
    for (uint i=0; i<4; ++i){
        this->coeff[i] /= x;  
    }
    return *this; 
}

Matrix2x2 abs(const Matrix2x2 & M) {
    Matrix2x2 temp;
    for (uint i=0; i<4; ++i){
        temp.coeff[i] = std::abs(M.coeff[i]);  
    }
    return temp; 
}
	
}//namespace
