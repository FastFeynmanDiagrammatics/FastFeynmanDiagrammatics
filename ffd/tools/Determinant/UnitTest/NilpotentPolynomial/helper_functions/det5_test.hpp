namespace ffd::determinant::unit_test{

#ifdef FFD_DETERMINANT_FIVE_UNIT_TEST_FLAG
  
  template<typename P_t>

  auto det5_test(P_t P){
    return std::vector<long double>{P[4][0]*P[8][0]*P[12][0]*P[16][0]*P[20][0] - 
   P[3][0]*P[9][0]*P[12][0]*P[16][0]*P[20][0] - 
   P[4][0]*P[7][0]*P[13][0]*P[16][0]*P[20][0] + 
   P[2][0]*P[9][0]*P[13][0]*P[16][0]*P[20][0] + 
   P[3][0]*P[7][0]*P[14][0]*P[16][0]*P[20][0] - 
   P[2][0]*P[8][0]*P[14][0]*P[16][0]*P[20][0] - 
   P[4][0]*P[8][0]*P[11][0]*P[17][0]*P[20][0] + 
   P[3][0]*P[9][0]*P[11][0]*P[17][0]*P[20][0] + 
   P[4][0]*P[6][0]*P[13][0]*P[17][0]*P[20][0] - 
   P[1][0]*P[9][0]*P[13][0]*P[17][0]*P[20][0] - 
   P[3][0]*P[6][0]*P[14][0]*P[17][0]*P[20][0] + 
   P[1][0]*P[8][0]*P[14][0]*P[17][0]*P[20][0] + 
   P[4][0]*P[7][0]*P[11][0]*P[18][0]*P[20][0] - 
   P[2][0]*P[9][0]*P[11][0]*P[18][0]*P[20][0] - 
   P[4][0]*P[6][0]*P[12][0]*P[18][0]*P[20][0] + 
   P[1][0]*P[9][0]*P[12][0]*P[18][0]*P[20][0] + 
   P[2][0]*P[6][0]*P[14][0]*P[18][0]*P[20][0] - 
   P[1][0]*P[7][0]*P[14][0]*P[18][0]*P[20][0] - 
   P[3][0]*P[7][0]*P[11][0]*P[19][0]*P[20][0] + 
   P[2][0]*P[8][0]*P[11][0]*P[19][0]*P[20][0] + 
   P[3][0]*P[6][0]*P[12][0]*P[19][0]*P[20][0] - 
   P[1][0]*P[8][0]*P[12][0]*P[19][0]*P[20][0] - 
   P[2][0]*P[6][0]*P[13][0]*P[19][0]*P[20][0] + 
   P[1][0]*P[7][0]*P[13][0]*P[19][0]*P[20][0] - 
   P[4][0]*P[8][0]*P[12][0]*P[15][0]*P[21][0] + 
   P[3][0]*P[9][0]*P[12][0]*P[15][0]*P[21][0] + 
   P[4][0]*P[7][0]*P[13][0]*P[15][0]*P[21][0] - 
   P[2][0]*P[9][0]*P[13][0]*P[15][0]*P[21][0] - 
   P[3][0]*P[7][0]*P[14][0]*P[15][0]*P[21][0] + 
   P[2][0]*P[8][0]*P[14][0]*P[15][0]*P[21][0] + 
   P[4][0]*P[8][0]*P[10][0]*P[17][0]*P[21][0] - 
   P[3][0]*P[9][0]*P[10][0]*P[17][0]*P[21][0] - 
   P[4][0]*P[5][0]*P[13][0]*P[17][0]*P[21][0] + 
   P[0][0]*P[9][0]*P[13][0]*P[17][0]*P[21][0] + 
   P[3][0]*P[5][0]*P[14][0]*P[17][0]*P[21][0] - 
   P[0][0]*P[8][0]*P[14][0]*P[17][0]*P[21][0] - 
   P[4][0]*P[7][0]*P[10][0]*P[18][0]*P[21][0] + 
   P[2][0]*P[9][0]*P[10][0]*P[18][0]*P[21][0] + 
   P[4][0]*P[5][0]*P[12][0]*P[18][0]*P[21][0] - 
   P[0][0]*P[9][0]*P[12][0]*P[18][0]*P[21][0] - 
   P[2][0]*P[5][0]*P[14][0]*P[18][0]*P[21][0] + 
   P[0][0]*P[7][0]*P[14][0]*P[18][0]*P[21][0] + 
   P[3][0]*P[7][0]*P[10][0]*P[19][0]*P[21][0] - 
   P[2][0]*P[8][0]*P[10][0]*P[19][0]*P[21][0] - 
   P[3][0]*P[5][0]*P[12][0]*P[19][0]*P[21][0] + 
   P[0][0]*P[8][0]*P[12][0]*P[19][0]*P[21][0] + 
   P[2][0]*P[5][0]*P[13][0]*P[19][0]*P[21][0] - 
   P[0][0]*P[7][0]*P[13][0]*P[19][0]*P[21][0] + 
   P[4][0]*P[8][0]*P[11][0]*P[15][0]*P[22][0] - 
   P[3][0]*P[9][0]*P[11][0]*P[15][0]*P[22][0] - 
   P[4][0]*P[6][0]*P[13][0]*P[15][0]*P[22][0] + 
   P[1][0]*P[9][0]*P[13][0]*P[15][0]*P[22][0] + 
   P[3][0]*P[6][0]*P[14][0]*P[15][0]*P[22][0] - 
   P[1][0]*P[8][0]*P[14][0]*P[15][0]*P[22][0] - 
   P[4][0]*P[8][0]*P[10][0]*P[16][0]*P[22][0] + 
   P[3][0]*P[9][0]*P[10][0]*P[16][0]*P[22][0] + 
   P[4][0]*P[5][0]*P[13][0]*P[16][0]*P[22][0] - 
   P[0][0]*P[9][0]*P[13][0]*P[16][0]*P[22][0] - 
   P[3][0]*P[5][0]*P[14][0]*P[16][0]*P[22][0] + 
   P[0][0]*P[8][0]*P[14][0]*P[16][0]*P[22][0] + 
   P[4][0]*P[6][0]*P[10][0]*P[18][0]*P[22][0] - 
   P[1][0]*P[9][0]*P[10][0]*P[18][0]*P[22][0] - 
   P[4][0]*P[5][0]*P[11][0]*P[18][0]*P[22][0] + 
   P[0][0]*P[9][0]*P[11][0]*P[18][0]*P[22][0] + 
   P[1][0]*P[5][0]*P[14][0]*P[18][0]*P[22][0] - 
   P[0][0]*P[6][0]*P[14][0]*P[18][0]*P[22][0] - 
   P[3][0]*P[6][0]*P[10][0]*P[19][0]*P[22][0] + 
   P[1][0]*P[8][0]*P[10][0]*P[19][0]*P[22][0] + 
   P[3][0]*P[5][0]*P[11][0]*P[19][0]*P[22][0] - 
   P[0][0]*P[8][0]*P[11][0]*P[19][0]*P[22][0] - 
   P[1][0]*P[5][0]*P[13][0]*P[19][0]*P[22][0] + 
   P[0][0]*P[6][0]*P[13][0]*P[19][0]*P[22][0] - 
   P[4][0]*P[7][0]*P[11][0]*P[15][0]*P[23][0] + 
   P[2][0]*P[9][0]*P[11][0]*P[15][0]*P[23][0] + 
   P[4][0]*P[6][0]*P[12][0]*P[15][0]*P[23][0] - 
   P[1][0]*P[9][0]*P[12][0]*P[15][0]*P[23][0] - 
   P[2][0]*P[6][0]*P[14][0]*P[15][0]*P[23][0] + 
   P[1][0]*P[7][0]*P[14][0]*P[15][0]*P[23][0] + 
   P[4][0]*P[7][0]*P[10][0]*P[16][0]*P[23][0] - 
   P[2][0]*P[9][0]*P[10][0]*P[16][0]*P[23][0] - 
   P[4][0]*P[5][0]*P[12][0]*P[16][0]*P[23][0] + 
   P[0][0]*P[9][0]*P[12][0]*P[16][0]*P[23][0] + 
   P[2][0]*P[5][0]*P[14][0]*P[16][0]*P[23][0] - 
   P[0][0]*P[7][0]*P[14][0]*P[16][0]*P[23][0] - 
   P[4][0]*P[6][0]*P[10][0]*P[17][0]*P[23][0] + 
   P[1][0]*P[9][0]*P[10][0]*P[17][0]*P[23][0] + 
   P[4][0]*P[5][0]*P[11][0]*P[17][0]*P[23][0] - 
   P[0][0]*P[9][0]*P[11][0]*P[17][0]*P[23][0] - 
   P[1][0]*P[5][0]*P[14][0]*P[17][0]*P[23][0] + 
   P[0][0]*P[6][0]*P[14][0]*P[17][0]*P[23][0] + 
   P[2][0]*P[6][0]*P[10][0]*P[19][0]*P[23][0] - 
   P[1][0]*P[7][0]*P[10][0]*P[19][0]*P[23][0] - 
   P[2][0]*P[5][0]*P[11][0]*P[19][0]*P[23][0] + 
   P[0][0]*P[7][0]*P[11][0]*P[19][0]*P[23][0] + 
   P[1][0]*P[5][0]*P[12][0]*P[19][0]*P[23][0] - 
   P[0][0]*P[6][0]*P[12][0]*P[19][0]*P[23][0] + 
   P[3][0]*P[7][0]*P[11][0]*P[15][0]*P[24][0] - 
   P[2][0]*P[8][0]*P[11][0]*P[15][0]*P[24][0] - 
   P[3][0]*P[6][0]*P[12][0]*P[15][0]*P[24][0] + 
   P[1][0]*P[8][0]*P[12][0]*P[15][0]*P[24][0] + 
   P[2][0]*P[6][0]*P[13][0]*P[15][0]*P[24][0] - 
   P[1][0]*P[7][0]*P[13][0]*P[15][0]*P[24][0] - 
   P[3][0]*P[7][0]*P[10][0]*P[16][0]*P[24][0] + 
   P[2][0]*P[8][0]*P[10][0]*P[16][0]*P[24][0] + 
   P[3][0]*P[5][0]*P[12][0]*P[16][0]*P[24][0] - 
   P[0][0]*P[8][0]*P[12][0]*P[16][0]*P[24][0] - 
   P[2][0]*P[5][0]*P[13][0]*P[16][0]*P[24][0] + 
   P[0][0]*P[7][0]*P[13][0]*P[16][0]*P[24][0] + 
   P[3][0]*P[6][0]*P[10][0]*P[17][0]*P[24][0] - 
   P[1][0]*P[8][0]*P[10][0]*P[17][0]*P[24][0] - 
   P[3][0]*P[5][0]*P[11][0]*P[17][0]*P[24][0] + 
   P[0][0]*P[8][0]*P[11][0]*P[17][0]*P[24][0] + 
   P[1][0]*P[5][0]*P[13][0]*P[17][0]*P[24][0] - 
   P[0][0]*P[6][0]*P[13][0]*P[17][0]*P[24][0] - 
   P[2][0]*P[6][0]*P[10][0]*P[18][0]*P[24][0] + 
   P[1][0]*P[7][0]*P[10][0]*P[18][0]*P[24][0] + 
   P[2][0]*P[5][0]*P[11][0]*P[18][0]*P[24][0] - 
   P[0][0]*P[7][0]*P[11][0]*P[18][0]*P[24][0] - 
   P[1][0]*P[5][0]*P[12][0]*P[18][0]*P[24][0] + 
				      P[0][0]*P[6][0]*P[12][0]*P[18][0]*P[24][0],
				        P[4][1]*P[8][0]*P[12][0]*P[16][0]*P[20][0] + 
   P[4][0]*P[8][1]*P[12][0]*P[16][0]*P[20][0] - 
   P[3][1]*P[9][0]*P[12][0]*P[16][0]*P[20][0] - 
   P[3][0]*P[9][1]*P[12][0]*P[16][0]*P[20][0] + 
   P[4][0]*P[8][0]*P[12][1]*P[16][0]*P[20][0] - 
   P[3][0]*P[9][0]*P[12][1]*P[16][0]*P[20][0] - 
   P[4][1]*P[7][0]*P[13][0]*P[16][0]*P[20][0] - 
   P[4][0]*P[7][1]*P[13][0]*P[16][0]*P[20][0] + 
   P[2][1]*P[9][0]*P[13][0]*P[16][0]*P[20][0] + 
   P[2][0]*P[9][1]*P[13][0]*P[16][0]*P[20][0] - 
   P[4][0]*P[7][0]*P[13][1]*P[16][0]*P[20][0] + 
   P[2][0]*P[9][0]*P[13][1]*P[16][0]*P[20][0] + 
   P[3][1]*P[7][0]*P[14][0]*P[16][0]*P[20][0] + 
   P[3][0]*P[7][1]*P[14][0]*P[16][0]*P[20][0] - 
   P[2][1]*P[8][0]*P[14][0]*P[16][0]*P[20][0] - 
   P[2][0]*P[8][1]*P[14][0]*P[16][0]*P[20][0] + 
   P[3][0]*P[7][0]*P[14][1]*P[16][0]*P[20][0] - 
   P[2][0]*P[8][0]*P[14][1]*P[16][0]*P[20][0] + 
   P[4][0]*P[8][0]*P[12][0]*P[16][1]*P[20][0] - 
   P[3][0]*P[9][0]*P[12][0]*P[16][1]*P[20][0] - 
   P[4][0]*P[7][0]*P[13][0]*P[16][1]*P[20][0] + 
   P[2][0]*P[9][0]*P[13][0]*P[16][1]*P[20][0] + 
   P[3][0]*P[7][0]*P[14][0]*P[16][1]*P[20][0] - 
   P[2][0]*P[8][0]*P[14][0]*P[16][1]*P[20][0] - 
   P[4][1]*P[8][0]*P[11][0]*P[17][0]*P[20][0] - 
   P[4][0]*P[8][1]*P[11][0]*P[17][0]*P[20][0] + 
   P[3][1]*P[9][0]*P[11][0]*P[17][0]*P[20][0] + 
   P[3][0]*P[9][1]*P[11][0]*P[17][0]*P[20][0] - 
   P[4][0]*P[8][0]*P[11][1]*P[17][0]*P[20][0] + 
   P[3][0]*P[9][0]*P[11][1]*P[17][0]*P[20][0] + 
   P[4][1]*P[6][0]*P[13][0]*P[17][0]*P[20][0] + 
   P[4][0]*P[6][1]*P[13][0]*P[17][0]*P[20][0] - 
   P[1][1]*P[9][0]*P[13][0]*P[17][0]*P[20][0] - 
   P[1][0]*P[9][1]*P[13][0]*P[17][0]*P[20][0] + 
   P[4][0]*P[6][0]*P[13][1]*P[17][0]*P[20][0] - 
   P[1][0]*P[9][0]*P[13][1]*P[17][0]*P[20][0] - 
   P[3][1]*P[6][0]*P[14][0]*P[17][0]*P[20][0] - 
   P[3][0]*P[6][1]*P[14][0]*P[17][0]*P[20][0] + 
   P[1][1]*P[8][0]*P[14][0]*P[17][0]*P[20][0] + 
   P[1][0]*P[8][1]*P[14][0]*P[17][0]*P[20][0] - 
   P[3][0]*P[6][0]*P[14][1]*P[17][0]*P[20][0] + 
   P[1][0]*P[8][0]*P[14][1]*P[17][0]*P[20][0] - 
   P[4][0]*P[8][0]*P[11][0]*P[17][1]*P[20][0] + 
   P[3][0]*P[9][0]*P[11][0]*P[17][1]*P[20][0] + 
   P[4][0]*P[6][0]*P[13][0]*P[17][1]*P[20][0] - 
   P[1][0]*P[9][0]*P[13][0]*P[17][1]*P[20][0] - 
   P[3][0]*P[6][0]*P[14][0]*P[17][1]*P[20][0] + 
   P[1][0]*P[8][0]*P[14][0]*P[17][1]*P[20][0] + 
   P[4][1]*P[7][0]*P[11][0]*P[18][0]*P[20][0] + 
   P[4][0]*P[7][1]*P[11][0]*P[18][0]*P[20][0] - 
   P[2][1]*P[9][0]*P[11][0]*P[18][0]*P[20][0] - 
   P[2][0]*P[9][1]*P[11][0]*P[18][0]*P[20][0] + 
   P[4][0]*P[7][0]*P[11][1]*P[18][0]*P[20][0] - 
   P[2][0]*P[9][0]*P[11][1]*P[18][0]*P[20][0] - 
   P[4][1]*P[6][0]*P[12][0]*P[18][0]*P[20][0] - 
   P[4][0]*P[6][1]*P[12][0]*P[18][0]*P[20][0] + 
   P[1][1]*P[9][0]*P[12][0]*P[18][0]*P[20][0] + 
   P[1][0]*P[9][1]*P[12][0]*P[18][0]*P[20][0] - 
   P[4][0]*P[6][0]*P[12][1]*P[18][0]*P[20][0] + 
   P[1][0]*P[9][0]*P[12][1]*P[18][0]*P[20][0] + 
   P[2][1]*P[6][0]*P[14][0]*P[18][0]*P[20][0] + 
   P[2][0]*P[6][1]*P[14][0]*P[18][0]*P[20][0] - 
   P[1][1]*P[7][0]*P[14][0]*P[18][0]*P[20][0] - 
   P[1][0]*P[7][1]*P[14][0]*P[18][0]*P[20][0] + 
   P[2][0]*P[6][0]*P[14][1]*P[18][0]*P[20][0] - 
   P[1][0]*P[7][0]*P[14][1]*P[18][0]*P[20][0] + 
   P[4][0]*P[7][0]*P[11][0]*P[18][1]*P[20][0] - 
   P[2][0]*P[9][0]*P[11][0]*P[18][1]*P[20][0] - 
   P[4][0]*P[6][0]*P[12][0]*P[18][1]*P[20][0] + 
   P[1][0]*P[9][0]*P[12][0]*P[18][1]*P[20][0] + 
   P[2][0]*P[6][0]*P[14][0]*P[18][1]*P[20][0] - 
   P[1][0]*P[7][0]*P[14][0]*P[18][1]*P[20][0] - 
   P[3][1]*P[7][0]*P[11][0]*P[19][0]*P[20][0] - 
   P[3][0]*P[7][1]*P[11][0]*P[19][0]*P[20][0] + 
   P[2][1]*P[8][0]*P[11][0]*P[19][0]*P[20][0] + 
   P[2][0]*P[8][1]*P[11][0]*P[19][0]*P[20][0] - 
   P[3][0]*P[7][0]*P[11][1]*P[19][0]*P[20][0] + 
   P[2][0]*P[8][0]*P[11][1]*P[19][0]*P[20][0] + 
   P[3][1]*P[6][0]*P[12][0]*P[19][0]*P[20][0] + 
   P[3][0]*P[6][1]*P[12][0]*P[19][0]*P[20][0] - 
   P[1][1]*P[8][0]*P[12][0]*P[19][0]*P[20][0] - 
   P[1][0]*P[8][1]*P[12][0]*P[19][0]*P[20][0] + 
   P[3][0]*P[6][0]*P[12][1]*P[19][0]*P[20][0] - 
   P[1][0]*P[8][0]*P[12][1]*P[19][0]*P[20][0] - 
   P[2][1]*P[6][0]*P[13][0]*P[19][0]*P[20][0] - 
   P[2][0]*P[6][1]*P[13][0]*P[19][0]*P[20][0] + 
   P[1][1]*P[7][0]*P[13][0]*P[19][0]*P[20][0] + 
   P[1][0]*P[7][1]*P[13][0]*P[19][0]*P[20][0] - 
   P[2][0]*P[6][0]*P[13][1]*P[19][0]*P[20][0] + 
   P[1][0]*P[7][0]*P[13][1]*P[19][0]*P[20][0] - 
   P[3][0]*P[7][0]*P[11][0]*P[19][1]*P[20][0] + 
   P[2][0]*P[8][0]*P[11][0]*P[19][1]*P[20][0] + 
   P[3][0]*P[6][0]*P[12][0]*P[19][1]*P[20][0] - 
   P[1][0]*P[8][0]*P[12][0]*P[19][1]*P[20][0] - 
   P[2][0]*P[6][0]*P[13][0]*P[19][1]*P[20][0] + 
   P[1][0]*P[7][0]*P[13][0]*P[19][1]*P[20][0] + 
   P[4][0]*P[8][0]*P[12][0]*P[16][0]*P[20][1] - 
   P[3][0]*P[9][0]*P[12][0]*P[16][0]*P[20][1] - 
   P[4][0]*P[7][0]*P[13][0]*P[16][0]*P[20][1] + 
   P[2][0]*P[9][0]*P[13][0]*P[16][0]*P[20][1] + 
   P[3][0]*P[7][0]*P[14][0]*P[16][0]*P[20][1] - 
   P[2][0]*P[8][0]*P[14][0]*P[16][0]*P[20][1] - 
   P[4][0]*P[8][0]*P[11][0]*P[17][0]*P[20][1] + 
   P[3][0]*P[9][0]*P[11][0]*P[17][0]*P[20][1] + 
   P[4][0]*P[6][0]*P[13][0]*P[17][0]*P[20][1] - 
   P[1][0]*P[9][0]*P[13][0]*P[17][0]*P[20][1] - 
   P[3][0]*P[6][0]*P[14][0]*P[17][0]*P[20][1] + 
   P[1][0]*P[8][0]*P[14][0]*P[17][0]*P[20][1] + 
   P[4][0]*P[7][0]*P[11][0]*P[18][0]*P[20][1] - 
   P[2][0]*P[9][0]*P[11][0]*P[18][0]*P[20][1] - 
   P[4][0]*P[6][0]*P[12][0]*P[18][0]*P[20][1] + 
   P[1][0]*P[9][0]*P[12][0]*P[18][0]*P[20][1] + 
   P[2][0]*P[6][0]*P[14][0]*P[18][0]*P[20][1] - 
   P[1][0]*P[7][0]*P[14][0]*P[18][0]*P[20][1] - 
   P[3][0]*P[7][0]*P[11][0]*P[19][0]*P[20][1] + 
   P[2][0]*P[8][0]*P[11][0]*P[19][0]*P[20][1] + 
   P[3][0]*P[6][0]*P[12][0]*P[19][0]*P[20][1] - 
   P[1][0]*P[8][0]*P[12][0]*P[19][0]*P[20][1] - 
   P[2][0]*P[6][0]*P[13][0]*P[19][0]*P[20][1] + 
   P[1][0]*P[7][0]*P[13][0]*P[19][0]*P[20][1] - 
   P[4][1]*P[8][0]*P[12][0]*P[15][0]*P[21][0] - 
   P[4][0]*P[8][1]*P[12][0]*P[15][0]*P[21][0] + 
   P[3][1]*P[9][0]*P[12][0]*P[15][0]*P[21][0] + 
   P[3][0]*P[9][1]*P[12][0]*P[15][0]*P[21][0] - 
   P[4][0]*P[8][0]*P[12][1]*P[15][0]*P[21][0] + 
   P[3][0]*P[9][0]*P[12][1]*P[15][0]*P[21][0] + 
   P[4][1]*P[7][0]*P[13][0]*P[15][0]*P[21][0] + 
   P[4][0]*P[7][1]*P[13][0]*P[15][0]*P[21][0] - 
   P[2][1]*P[9][0]*P[13][0]*P[15][0]*P[21][0] - 
   P[2][0]*P[9][1]*P[13][0]*P[15][0]*P[21][0] + 
   P[4][0]*P[7][0]*P[13][1]*P[15][0]*P[21][0] - 
   P[2][0]*P[9][0]*P[13][1]*P[15][0]*P[21][0] - 
   P[3][1]*P[7][0]*P[14][0]*P[15][0]*P[21][0] - 
   P[3][0]*P[7][1]*P[14][0]*P[15][0]*P[21][0] + 
   P[2][1]*P[8][0]*P[14][0]*P[15][0]*P[21][0] + 
   P[2][0]*P[8][1]*P[14][0]*P[15][0]*P[21][0] - 
   P[3][0]*P[7][0]*P[14][1]*P[15][0]*P[21][0] + 
   P[2][0]*P[8][0]*P[14][1]*P[15][0]*P[21][0] - 
   P[4][0]*P[8][0]*P[12][0]*P[15][1]*P[21][0] + 
   P[3][0]*P[9][0]*P[12][0]*P[15][1]*P[21][0] + 
   P[4][0]*P[7][0]*P[13][0]*P[15][1]*P[21][0] - 
   P[2][0]*P[9][0]*P[13][0]*P[15][1]*P[21][0] - 
   P[3][0]*P[7][0]*P[14][0]*P[15][1]*P[21][0] + 
   P[2][0]*P[8][0]*P[14][0]*P[15][1]*P[21][0] + 
   P[4][1]*P[8][0]*P[10][0]*P[17][0]*P[21][0] + 
   P[4][0]*P[8][1]*P[10][0]*P[17][0]*P[21][0] - 
   P[3][1]*P[9][0]*P[10][0]*P[17][0]*P[21][0] - 
   P[3][0]*P[9][1]*P[10][0]*P[17][0]*P[21][0] + 
   P[4][0]*P[8][0]*P[10][1]*P[17][0]*P[21][0] - 
   P[3][0]*P[9][0]*P[10][1]*P[17][0]*P[21][0] - 
   P[4][1]*P[5][0]*P[13][0]*P[17][0]*P[21][0] - 
   P[4][0]*P[5][1]*P[13][0]*P[17][0]*P[21][0] + 
   P[0][1]*P[9][0]*P[13][0]*P[17][0]*P[21][0] + 
   P[0][0]*P[9][1]*P[13][0]*P[17][0]*P[21][0] - 
   P[4][0]*P[5][0]*P[13][1]*P[17][0]*P[21][0] + 
   P[0][0]*P[9][0]*P[13][1]*P[17][0]*P[21][0] + 
   P[3][1]*P[5][0]*P[14][0]*P[17][0]*P[21][0] + 
   P[3][0]*P[5][1]*P[14][0]*P[17][0]*P[21][0] - 
   P[0][1]*P[8][0]*P[14][0]*P[17][0]*P[21][0] - 
   P[0][0]*P[8][1]*P[14][0]*P[17][0]*P[21][0] + 
   P[3][0]*P[5][0]*P[14][1]*P[17][0]*P[21][0] - 
   P[0][0]*P[8][0]*P[14][1]*P[17][0]*P[21][0] + 
   P[4][0]*P[8][0]*P[10][0]*P[17][1]*P[21][0] - 
   P[3][0]*P[9][0]*P[10][0]*P[17][1]*P[21][0] - 
   P[4][0]*P[5][0]*P[13][0]*P[17][1]*P[21][0] + 
   P[0][0]*P[9][0]*P[13][0]*P[17][1]*P[21][0] + 
   P[3][0]*P[5][0]*P[14][0]*P[17][1]*P[21][0] - 
   P[0][0]*P[8][0]*P[14][0]*P[17][1]*P[21][0] - 
   P[4][1]*P[7][0]*P[10][0]*P[18][0]*P[21][0] - 
   P[4][0]*P[7][1]*P[10][0]*P[18][0]*P[21][0] + 
   P[2][1]*P[9][0]*P[10][0]*P[18][0]*P[21][0] + 
   P[2][0]*P[9][1]*P[10][0]*P[18][0]*P[21][0] - 
   P[4][0]*P[7][0]*P[10][1]*P[18][0]*P[21][0] + 
   P[2][0]*P[9][0]*P[10][1]*P[18][0]*P[21][0] + 
   P[4][1]*P[5][0]*P[12][0]*P[18][0]*P[21][0] + 
   P[4][0]*P[5][1]*P[12][0]*P[18][0]*P[21][0] - 
   P[0][1]*P[9][0]*P[12][0]*P[18][0]*P[21][0] - 
   P[0][0]*P[9][1]*P[12][0]*P[18][0]*P[21][0] + 
   P[4][0]*P[5][0]*P[12][1]*P[18][0]*P[21][0] - 
   P[0][0]*P[9][0]*P[12][1]*P[18][0]*P[21][0] - 
   P[2][1]*P[5][0]*P[14][0]*P[18][0]*P[21][0] - 
   P[2][0]*P[5][1]*P[14][0]*P[18][0]*P[21][0] + 
   P[0][1]*P[7][0]*P[14][0]*P[18][0]*P[21][0] + 
   P[0][0]*P[7][1]*P[14][0]*P[18][0]*P[21][0] - 
   P[2][0]*P[5][0]*P[14][1]*P[18][0]*P[21][0] + 
   P[0][0]*P[7][0]*P[14][1]*P[18][0]*P[21][0] - 
   P[4][0]*P[7][0]*P[10][0]*P[18][1]*P[21][0] + 
   P[2][0]*P[9][0]*P[10][0]*P[18][1]*P[21][0] + 
   P[4][0]*P[5][0]*P[12][0]*P[18][1]*P[21][0] - 
   P[0][0]*P[9][0]*P[12][0]*P[18][1]*P[21][0] - 
   P[2][0]*P[5][0]*P[14][0]*P[18][1]*P[21][0] + 
   P[0][0]*P[7][0]*P[14][0]*P[18][1]*P[21][0] + 
   P[3][1]*P[7][0]*P[10][0]*P[19][0]*P[21][0] + 
   P[3][0]*P[7][1]*P[10][0]*P[19][0]*P[21][0] - 
   P[2][1]*P[8][0]*P[10][0]*P[19][0]*P[21][0] - 
   P[2][0]*P[8][1]*P[10][0]*P[19][0]*P[21][0] + 
   P[3][0]*P[7][0]*P[10][1]*P[19][0]*P[21][0] - 
   P[2][0]*P[8][0]*P[10][1]*P[19][0]*P[21][0] - 
   P[3][1]*P[5][0]*P[12][0]*P[19][0]*P[21][0] - 
   P[3][0]*P[5][1]*P[12][0]*P[19][0]*P[21][0] + 
   P[0][1]*P[8][0]*P[12][0]*P[19][0]*P[21][0] + 
   P[0][0]*P[8][1]*P[12][0]*P[19][0]*P[21][0] - 
   P[3][0]*P[5][0]*P[12][1]*P[19][0]*P[21][0] + 
   P[0][0]*P[8][0]*P[12][1]*P[19][0]*P[21][0] + 
   P[2][1]*P[5][0]*P[13][0]*P[19][0]*P[21][0] + 
   P[2][0]*P[5][1]*P[13][0]*P[19][0]*P[21][0] - 
   P[0][1]*P[7][0]*P[13][0]*P[19][0]*P[21][0] - 
   P[0][0]*P[7][1]*P[13][0]*P[19][0]*P[21][0] + 
   P[2][0]*P[5][0]*P[13][1]*P[19][0]*P[21][0] - 
   P[0][0]*P[7][0]*P[13][1]*P[19][0]*P[21][0] + 
   P[3][0]*P[7][0]*P[10][0]*P[19][1]*P[21][0] - 
   P[2][0]*P[8][0]*P[10][0]*P[19][1]*P[21][0] - 
   P[3][0]*P[5][0]*P[12][0]*P[19][1]*P[21][0] + 
   P[0][0]*P[8][0]*P[12][0]*P[19][1]*P[21][0] + 
   P[2][0]*P[5][0]*P[13][0]*P[19][1]*P[21][0] - 
   P[0][0]*P[7][0]*P[13][0]*P[19][1]*P[21][0] - 
   P[4][0]*P[8][0]*P[12][0]*P[15][0]*P[21][1] + 
   P[3][0]*P[9][0]*P[12][0]*P[15][0]*P[21][1] + 
   P[4][0]*P[7][0]*P[13][0]*P[15][0]*P[21][1] - 
   P[2][0]*P[9][0]*P[13][0]*P[15][0]*P[21][1] - 
   P[3][0]*P[7][0]*P[14][0]*P[15][0]*P[21][1] + 
   P[2][0]*P[8][0]*P[14][0]*P[15][0]*P[21][1] + 
   P[4][0]*P[8][0]*P[10][0]*P[17][0]*P[21][1] - 
   P[3][0]*P[9][0]*P[10][0]*P[17][0]*P[21][1] - 
   P[4][0]*P[5][0]*P[13][0]*P[17][0]*P[21][1] + 
   P[0][0]*P[9][0]*P[13][0]*P[17][0]*P[21][1] + 
   P[3][0]*P[5][0]*P[14][0]*P[17][0]*P[21][1] - 
   P[0][0]*P[8][0]*P[14][0]*P[17][0]*P[21][1] - 
   P[4][0]*P[7][0]*P[10][0]*P[18][0]*P[21][1] + 
   P[2][0]*P[9][0]*P[10][0]*P[18][0]*P[21][1] + 
   P[4][0]*P[5][0]*P[12][0]*P[18][0]*P[21][1] - 
   P[0][0]*P[9][0]*P[12][0]*P[18][0]*P[21][1] - 
   P[2][0]*P[5][0]*P[14][0]*P[18][0]*P[21][1] + 
   P[0][0]*P[7][0]*P[14][0]*P[18][0]*P[21][1] + 
   P[3][0]*P[7][0]*P[10][0]*P[19][0]*P[21][1] - 
   P[2][0]*P[8][0]*P[10][0]*P[19][0]*P[21][1] - 
   P[3][0]*P[5][0]*P[12][0]*P[19][0]*P[21][1] + 
   P[0][0]*P[8][0]*P[12][0]*P[19][0]*P[21][1] + 
   P[2][0]*P[5][0]*P[13][0]*P[19][0]*P[21][1] - 
   P[0][0]*P[7][0]*P[13][0]*P[19][0]*P[21][1] + 
   P[4][1]*P[8][0]*P[11][0]*P[15][0]*P[22][0] + 
   P[4][0]*P[8][1]*P[11][0]*P[15][0]*P[22][0] - 
   P[3][1]*P[9][0]*P[11][0]*P[15][0]*P[22][0] - 
   P[3][0]*P[9][1]*P[11][0]*P[15][0]*P[22][0] + 
   P[4][0]*P[8][0]*P[11][1]*P[15][0]*P[22][0] - 
   P[3][0]*P[9][0]*P[11][1]*P[15][0]*P[22][0] - 
   P[4][1]*P[6][0]*P[13][0]*P[15][0]*P[22][0] - 
   P[4][0]*P[6][1]*P[13][0]*P[15][0]*P[22][0] + 
   P[1][1]*P[9][0]*P[13][0]*P[15][0]*P[22][0] + 
   P[1][0]*P[9][1]*P[13][0]*P[15][0]*P[22][0] - 
   P[4][0]*P[6][0]*P[13][1]*P[15][0]*P[22][0] + 
   P[1][0]*P[9][0]*P[13][1]*P[15][0]*P[22][0] + 
   P[3][1]*P[6][0]*P[14][0]*P[15][0]*P[22][0] + 
   P[3][0]*P[6][1]*P[14][0]*P[15][0]*P[22][0] - 
   P[1][1]*P[8][0]*P[14][0]*P[15][0]*P[22][0] - 
   P[1][0]*P[8][1]*P[14][0]*P[15][0]*P[22][0] + 
   P[3][0]*P[6][0]*P[14][1]*P[15][0]*P[22][0] - 
   P[1][0]*P[8][0]*P[14][1]*P[15][0]*P[22][0] + 
   P[4][0]*P[8][0]*P[11][0]*P[15][1]*P[22][0] - 
   P[3][0]*P[9][0]*P[11][0]*P[15][1]*P[22][0] - 
   P[4][0]*P[6][0]*P[13][0]*P[15][1]*P[22][0] + 
   P[1][0]*P[9][0]*P[13][0]*P[15][1]*P[22][0] + 
   P[3][0]*P[6][0]*P[14][0]*P[15][1]*P[22][0] - 
   P[1][0]*P[8][0]*P[14][0]*P[15][1]*P[22][0] - 
   P[4][1]*P[8][0]*P[10][0]*P[16][0]*P[22][0] - 
   P[4][0]*P[8][1]*P[10][0]*P[16][0]*P[22][0] + 
   P[3][1]*P[9][0]*P[10][0]*P[16][0]*P[22][0] + 
   P[3][0]*P[9][1]*P[10][0]*P[16][0]*P[22][0] - 
   P[4][0]*P[8][0]*P[10][1]*P[16][0]*P[22][0] + 
   P[3][0]*P[9][0]*P[10][1]*P[16][0]*P[22][0] + 
   P[4][1]*P[5][0]*P[13][0]*P[16][0]*P[22][0] + 
   P[4][0]*P[5][1]*P[13][0]*P[16][0]*P[22][0] - 
   P[0][1]*P[9][0]*P[13][0]*P[16][0]*P[22][0] - 
   P[0][0]*P[9][1]*P[13][0]*P[16][0]*P[22][0] + 
   P[4][0]*P[5][0]*P[13][1]*P[16][0]*P[22][0] - 
   P[0][0]*P[9][0]*P[13][1]*P[16][0]*P[22][0] - 
   P[3][1]*P[5][0]*P[14][0]*P[16][0]*P[22][0] - 
   P[3][0]*P[5][1]*P[14][0]*P[16][0]*P[22][0] + 
   P[0][1]*P[8][0]*P[14][0]*P[16][0]*P[22][0] + 
   P[0][0]*P[8][1]*P[14][0]*P[16][0]*P[22][0] - 
   P[3][0]*P[5][0]*P[14][1]*P[16][0]*P[22][0] + 
   P[0][0]*P[8][0]*P[14][1]*P[16][0]*P[22][0] - 
   P[4][0]*P[8][0]*P[10][0]*P[16][1]*P[22][0] + 
   P[3][0]*P[9][0]*P[10][0]*P[16][1]*P[22][0] + 
   P[4][0]*P[5][0]*P[13][0]*P[16][1]*P[22][0] - 
   P[0][0]*P[9][0]*P[13][0]*P[16][1]*P[22][0] - 
   P[3][0]*P[5][0]*P[14][0]*P[16][1]*P[22][0] + 
   P[0][0]*P[8][0]*P[14][0]*P[16][1]*P[22][0] + 
   P[4][1]*P[6][0]*P[10][0]*P[18][0]*P[22][0] + 
   P[4][0]*P[6][1]*P[10][0]*P[18][0]*P[22][0] - 
   P[1][1]*P[9][0]*P[10][0]*P[18][0]*P[22][0] - 
   P[1][0]*P[9][1]*P[10][0]*P[18][0]*P[22][0] + 
   P[4][0]*P[6][0]*P[10][1]*P[18][0]*P[22][0] - 
   P[1][0]*P[9][0]*P[10][1]*P[18][0]*P[22][0] - 
   P[4][1]*P[5][0]*P[11][0]*P[18][0]*P[22][0] - 
   P[4][0]*P[5][1]*P[11][0]*P[18][0]*P[22][0] + 
   P[0][1]*P[9][0]*P[11][0]*P[18][0]*P[22][0] + 
   P[0][0]*P[9][1]*P[11][0]*P[18][0]*P[22][0] - 
   P[4][0]*P[5][0]*P[11][1]*P[18][0]*P[22][0] + 
   P[0][0]*P[9][0]*P[11][1]*P[18][0]*P[22][0] + 
   P[1][1]*P[5][0]*P[14][0]*P[18][0]*P[22][0] + 
   P[1][0]*P[5][1]*P[14][0]*P[18][0]*P[22][0] - 
   P[0][1]*P[6][0]*P[14][0]*P[18][0]*P[22][0] - 
   P[0][0]*P[6][1]*P[14][0]*P[18][0]*P[22][0] + 
   P[1][0]*P[5][0]*P[14][1]*P[18][0]*P[22][0] - 
   P[0][0]*P[6][0]*P[14][1]*P[18][0]*P[22][0] + 
   P[4][0]*P[6][0]*P[10][0]*P[18][1]*P[22][0] - 
   P[1][0]*P[9][0]*P[10][0]*P[18][1]*P[22][0] - 
   P[4][0]*P[5][0]*P[11][0]*P[18][1]*P[22][0] + 
   P[0][0]*P[9][0]*P[11][0]*P[18][1]*P[22][0] + 
   P[1][0]*P[5][0]*P[14][0]*P[18][1]*P[22][0] - 
   P[0][0]*P[6][0]*P[14][0]*P[18][1]*P[22][0] - 
   P[3][1]*P[6][0]*P[10][0]*P[19][0]*P[22][0] - 
   P[3][0]*P[6][1]*P[10][0]*P[19][0]*P[22][0] + 
   P[1][1]*P[8][0]*P[10][0]*P[19][0]*P[22][0] + 
   P[1][0]*P[8][1]*P[10][0]*P[19][0]*P[22][0] - 
   P[3][0]*P[6][0]*P[10][1]*P[19][0]*P[22][0] + 
   P[1][0]*P[8][0]*P[10][1]*P[19][0]*P[22][0] + 
   P[3][1]*P[5][0]*P[11][0]*P[19][0]*P[22][0] + 
   P[3][0]*P[5][1]*P[11][0]*P[19][0]*P[22][0] - 
   P[0][1]*P[8][0]*P[11][0]*P[19][0]*P[22][0] - 
   P[0][0]*P[8][1]*P[11][0]*P[19][0]*P[22][0] + 
   P[3][0]*P[5][0]*P[11][1]*P[19][0]*P[22][0] - 
   P[0][0]*P[8][0]*P[11][1]*P[19][0]*P[22][0] - 
   P[1][1]*P[5][0]*P[13][0]*P[19][0]*P[22][0] - 
   P[1][0]*P[5][1]*P[13][0]*P[19][0]*P[22][0] + 
   P[0][1]*P[6][0]*P[13][0]*P[19][0]*P[22][0] + 
   P[0][0]*P[6][1]*P[13][0]*P[19][0]*P[22][0] - 
   P[1][0]*P[5][0]*P[13][1]*P[19][0]*P[22][0] + 
   P[0][0]*P[6][0]*P[13][1]*P[19][0]*P[22][0] - 
   P[3][0]*P[6][0]*P[10][0]*P[19][1]*P[22][0] + 
   P[1][0]*P[8][0]*P[10][0]*P[19][1]*P[22][0] + 
   P[3][0]*P[5][0]*P[11][0]*P[19][1]*P[22][0] - 
   P[0][0]*P[8][0]*P[11][0]*P[19][1]*P[22][0] - 
   P[1][0]*P[5][0]*P[13][0]*P[19][1]*P[22][0] + 
   P[0][0]*P[6][0]*P[13][0]*P[19][1]*P[22][0] + 
   P[4][0]*P[8][0]*P[11][0]*P[15][0]*P[22][1] - 
   P[3][0]*P[9][0]*P[11][0]*P[15][0]*P[22][1] - 
   P[4][0]*P[6][0]*P[13][0]*P[15][0]*P[22][1] + 
   P[1][0]*P[9][0]*P[13][0]*P[15][0]*P[22][1] + 
   P[3][0]*P[6][0]*P[14][0]*P[15][0]*P[22][1] - 
   P[1][0]*P[8][0]*P[14][0]*P[15][0]*P[22][1] - 
   P[4][0]*P[8][0]*P[10][0]*P[16][0]*P[22][1] + 
   P[3][0]*P[9][0]*P[10][0]*P[16][0]*P[22][1] + 
   P[4][0]*P[5][0]*P[13][0]*P[16][0]*P[22][1] - 
   P[0][0]*P[9][0]*P[13][0]*P[16][0]*P[22][1] - 
   P[3][0]*P[5][0]*P[14][0]*P[16][0]*P[22][1] + 
   P[0][0]*P[8][0]*P[14][0]*P[16][0]*P[22][1] + 
   P[4][0]*P[6][0]*P[10][0]*P[18][0]*P[22][1] - 
   P[1][0]*P[9][0]*P[10][0]*P[18][0]*P[22][1] - 
   P[4][0]*P[5][0]*P[11][0]*P[18][0]*P[22][1] + 
   P[0][0]*P[9][0]*P[11][0]*P[18][0]*P[22][1] + 
   P[1][0]*P[5][0]*P[14][0]*P[18][0]*P[22][1] - 
   P[0][0]*P[6][0]*P[14][0]*P[18][0]*P[22][1] - 
   P[3][0]*P[6][0]*P[10][0]*P[19][0]*P[22][1] + 
   P[1][0]*P[8][0]*P[10][0]*P[19][0]*P[22][1] + 
   P[3][0]*P[5][0]*P[11][0]*P[19][0]*P[22][1] - 
   P[0][0]*P[8][0]*P[11][0]*P[19][0]*P[22][1] - 
   P[1][0]*P[5][0]*P[13][0]*P[19][0]*P[22][1] + 
   P[0][0]*P[6][0]*P[13][0]*P[19][0]*P[22][1] - 
   P[4][1]*P[7][0]*P[11][0]*P[15][0]*P[23][0] - 
   P[4][0]*P[7][1]*P[11][0]*P[15][0]*P[23][0] + 
   P[2][1]*P[9][0]*P[11][0]*P[15][0]*P[23][0] + 
   P[2][0]*P[9][1]*P[11][0]*P[15][0]*P[23][0] - 
   P[4][0]*P[7][0]*P[11][1]*P[15][0]*P[23][0] + 
   P[2][0]*P[9][0]*P[11][1]*P[15][0]*P[23][0] + 
   P[4][1]*P[6][0]*P[12][0]*P[15][0]*P[23][0] + 
   P[4][0]*P[6][1]*P[12][0]*P[15][0]*P[23][0] - 
   P[1][1]*P[9][0]*P[12][0]*P[15][0]*P[23][0] - 
   P[1][0]*P[9][1]*P[12][0]*P[15][0]*P[23][0] + 
   P[4][0]*P[6][0]*P[12][1]*P[15][0]*P[23][0] - 
   P[1][0]*P[9][0]*P[12][1]*P[15][0]*P[23][0] - 
   P[2][1]*P[6][0]*P[14][0]*P[15][0]*P[23][0] - 
   P[2][0]*P[6][1]*P[14][0]*P[15][0]*P[23][0] + 
   P[1][1]*P[7][0]*P[14][0]*P[15][0]*P[23][0] + 
   P[1][0]*P[7][1]*P[14][0]*P[15][0]*P[23][0] - 
   P[2][0]*P[6][0]*P[14][1]*P[15][0]*P[23][0] + 
   P[1][0]*P[7][0]*P[14][1]*P[15][0]*P[23][0] - 
   P[4][0]*P[7][0]*P[11][0]*P[15][1]*P[23][0] + 
   P[2][0]*P[9][0]*P[11][0]*P[15][1]*P[23][0] + 
   P[4][0]*P[6][0]*P[12][0]*P[15][1]*P[23][0] - 
   P[1][0]*P[9][0]*P[12][0]*P[15][1]*P[23][0] - 
   P[2][0]*P[6][0]*P[14][0]*P[15][1]*P[23][0] + 
   P[1][0]*P[7][0]*P[14][0]*P[15][1]*P[23][0] + 
   P[4][1]*P[7][0]*P[10][0]*P[16][0]*P[23][0] + 
   P[4][0]*P[7][1]*P[10][0]*P[16][0]*P[23][0] - 
   P[2][1]*P[9][0]*P[10][0]*P[16][0]*P[23][0] - 
   P[2][0]*P[9][1]*P[10][0]*P[16][0]*P[23][0] + 
   P[4][0]*P[7][0]*P[10][1]*P[16][0]*P[23][0] - 
   P[2][0]*P[9][0]*P[10][1]*P[16][0]*P[23][0] - 
   P[4][1]*P[5][0]*P[12][0]*P[16][0]*P[23][0] - 
   P[4][0]*P[5][1]*P[12][0]*P[16][0]*P[23][0] + 
   P[0][1]*P[9][0]*P[12][0]*P[16][0]*P[23][0] + 
   P[0][0]*P[9][1]*P[12][0]*P[16][0]*P[23][0] - 
   P[4][0]*P[5][0]*P[12][1]*P[16][0]*P[23][0] + 
   P[0][0]*P[9][0]*P[12][1]*P[16][0]*P[23][0] + 
   P[2][1]*P[5][0]*P[14][0]*P[16][0]*P[23][0] + 
   P[2][0]*P[5][1]*P[14][0]*P[16][0]*P[23][0] - 
   P[0][1]*P[7][0]*P[14][0]*P[16][0]*P[23][0] - 
   P[0][0]*P[7][1]*P[14][0]*P[16][0]*P[23][0] + 
   P[2][0]*P[5][0]*P[14][1]*P[16][0]*P[23][0] - 
   P[0][0]*P[7][0]*P[14][1]*P[16][0]*P[23][0] + 
   P[4][0]*P[7][0]*P[10][0]*P[16][1]*P[23][0] - 
   P[2][0]*P[9][0]*P[10][0]*P[16][1]*P[23][0] - 
   P[4][0]*P[5][0]*P[12][0]*P[16][1]*P[23][0] + 
   P[0][0]*P[9][0]*P[12][0]*P[16][1]*P[23][0] + 
   P[2][0]*P[5][0]*P[14][0]*P[16][1]*P[23][0] - 
   P[0][0]*P[7][0]*P[14][0]*P[16][1]*P[23][0] - 
   P[4][1]*P[6][0]*P[10][0]*P[17][0]*P[23][0] - 
   P[4][0]*P[6][1]*P[10][0]*P[17][0]*P[23][0] + 
   P[1][1]*P[9][0]*P[10][0]*P[17][0]*P[23][0] + 
   P[1][0]*P[9][1]*P[10][0]*P[17][0]*P[23][0] - 
   P[4][0]*P[6][0]*P[10][1]*P[17][0]*P[23][0] + 
   P[1][0]*P[9][0]*P[10][1]*P[17][0]*P[23][0] + 
   P[4][1]*P[5][0]*P[11][0]*P[17][0]*P[23][0] + 
   P[4][0]*P[5][1]*P[11][0]*P[17][0]*P[23][0] - 
   P[0][1]*P[9][0]*P[11][0]*P[17][0]*P[23][0] - 
   P[0][0]*P[9][1]*P[11][0]*P[17][0]*P[23][0] + 
   P[4][0]*P[5][0]*P[11][1]*P[17][0]*P[23][0] - 
   P[0][0]*P[9][0]*P[11][1]*P[17][0]*P[23][0] - 
   P[1][1]*P[5][0]*P[14][0]*P[17][0]*P[23][0] - 
   P[1][0]*P[5][1]*P[14][0]*P[17][0]*P[23][0] + 
   P[0][1]*P[6][0]*P[14][0]*P[17][0]*P[23][0] + 
   P[0][0]*P[6][1]*P[14][0]*P[17][0]*P[23][0] - 
   P[1][0]*P[5][0]*P[14][1]*P[17][0]*P[23][0] + 
   P[0][0]*P[6][0]*P[14][1]*P[17][0]*P[23][0] - 
   P[4][0]*P[6][0]*P[10][0]*P[17][1]*P[23][0] + 
   P[1][0]*P[9][0]*P[10][0]*P[17][1]*P[23][0] + 
   P[4][0]*P[5][0]*P[11][0]*P[17][1]*P[23][0] - 
   P[0][0]*P[9][0]*P[11][0]*P[17][1]*P[23][0] - 
   P[1][0]*P[5][0]*P[14][0]*P[17][1]*P[23][0] + 
   P[0][0]*P[6][0]*P[14][0]*P[17][1]*P[23][0] + 
   P[2][1]*P[6][0]*P[10][0]*P[19][0]*P[23][0] + 
   P[2][0]*P[6][1]*P[10][0]*P[19][0]*P[23][0] - 
   P[1][1]*P[7][0]*P[10][0]*P[19][0]*P[23][0] - 
   P[1][0]*P[7][1]*P[10][0]*P[19][0]*P[23][0] + 
   P[2][0]*P[6][0]*P[10][1]*P[19][0]*P[23][0] - 
   P[1][0]*P[7][0]*P[10][1]*P[19][0]*P[23][0] - 
   P[2][1]*P[5][0]*P[11][0]*P[19][0]*P[23][0] - 
   P[2][0]*P[5][1]*P[11][0]*P[19][0]*P[23][0] + 
   P[0][1]*P[7][0]*P[11][0]*P[19][0]*P[23][0] + 
   P[0][0]*P[7][1]*P[11][0]*P[19][0]*P[23][0] - 
   P[2][0]*P[5][0]*P[11][1]*P[19][0]*P[23][0] + 
   P[0][0]*P[7][0]*P[11][1]*P[19][0]*P[23][0] + 
   P[1][1]*P[5][0]*P[12][0]*P[19][0]*P[23][0] + 
   P[1][0]*P[5][1]*P[12][0]*P[19][0]*P[23][0] - 
   P[0][1]*P[6][0]*P[12][0]*P[19][0]*P[23][0] - 
   P[0][0]*P[6][1]*P[12][0]*P[19][0]*P[23][0] + 
   P[1][0]*P[5][0]*P[12][1]*P[19][0]*P[23][0] - 
   P[0][0]*P[6][0]*P[12][1]*P[19][0]*P[23][0] + 
   P[2][0]*P[6][0]*P[10][0]*P[19][1]*P[23][0] - 
   P[1][0]*P[7][0]*P[10][0]*P[19][1]*P[23][0] - 
   P[2][0]*P[5][0]*P[11][0]*P[19][1]*P[23][0] + 
   P[0][0]*P[7][0]*P[11][0]*P[19][1]*P[23][0] + 
   P[1][0]*P[5][0]*P[12][0]*P[19][1]*P[23][0] - 
   P[0][0]*P[6][0]*P[12][0]*P[19][1]*P[23][0] - 
   P[4][0]*P[7][0]*P[11][0]*P[15][0]*P[23][1] + 
   P[2][0]*P[9][0]*P[11][0]*P[15][0]*P[23][1] + 
   P[4][0]*P[6][0]*P[12][0]*P[15][0]*P[23][1] - 
   P[1][0]*P[9][0]*P[12][0]*P[15][0]*P[23][1] - 
   P[2][0]*P[6][0]*P[14][0]*P[15][0]*P[23][1] + 
   P[1][0]*P[7][0]*P[14][0]*P[15][0]*P[23][1] + 
   P[4][0]*P[7][0]*P[10][0]*P[16][0]*P[23][1] - 
   P[2][0]*P[9][0]*P[10][0]*P[16][0]*P[23][1] - 
   P[4][0]*P[5][0]*P[12][0]*P[16][0]*P[23][1] + 
   P[0][0]*P[9][0]*P[12][0]*P[16][0]*P[23][1] + 
   P[2][0]*P[5][0]*P[14][0]*P[16][0]*P[23][1] - 
   P[0][0]*P[7][0]*P[14][0]*P[16][0]*P[23][1] - 
   P[4][0]*P[6][0]*P[10][0]*P[17][0]*P[23][1] + 
   P[1][0]*P[9][0]*P[10][0]*P[17][0]*P[23][1] + 
   P[4][0]*P[5][0]*P[11][0]*P[17][0]*P[23][1] - 
   P[0][0]*P[9][0]*P[11][0]*P[17][0]*P[23][1] - 
   P[1][0]*P[5][0]*P[14][0]*P[17][0]*P[23][1] + 
   P[0][0]*P[6][0]*P[14][0]*P[17][0]*P[23][1] + 
   P[2][0]*P[6][0]*P[10][0]*P[19][0]*P[23][1] - 
   P[1][0]*P[7][0]*P[10][0]*P[19][0]*P[23][1] - 
   P[2][0]*P[5][0]*P[11][0]*P[19][0]*P[23][1] + 
   P[0][0]*P[7][0]*P[11][0]*P[19][0]*P[23][1] + 
   P[1][0]*P[5][0]*P[12][0]*P[19][0]*P[23][1] - 
   P[0][0]*P[6][0]*P[12][0]*P[19][0]*P[23][1] + 
   P[3][1]*P[7][0]*P[11][0]*P[15][0]*P[24][0] + 
   P[3][0]*P[7][1]*P[11][0]*P[15][0]*P[24][0] - 
   P[2][1]*P[8][0]*P[11][0]*P[15][0]*P[24][0] - 
   P[2][0]*P[8][1]*P[11][0]*P[15][0]*P[24][0] + 
   P[3][0]*P[7][0]*P[11][1]*P[15][0]*P[24][0] - 
   P[2][0]*P[8][0]*P[11][1]*P[15][0]*P[24][0] - 
   P[3][1]*P[6][0]*P[12][0]*P[15][0]*P[24][0] - 
   P[3][0]*P[6][1]*P[12][0]*P[15][0]*P[24][0] + 
   P[1][1]*P[8][0]*P[12][0]*P[15][0]*P[24][0] + 
   P[1][0]*P[8][1]*P[12][0]*P[15][0]*P[24][0] - 
   P[3][0]*P[6][0]*P[12][1]*P[15][0]*P[24][0] + 
   P[1][0]*P[8][0]*P[12][1]*P[15][0]*P[24][0] + 
   P[2][1]*P[6][0]*P[13][0]*P[15][0]*P[24][0] + 
   P[2][0]*P[6][1]*P[13][0]*P[15][0]*P[24][0] - 
   P[1][1]*P[7][0]*P[13][0]*P[15][0]*P[24][0] - 
   P[1][0]*P[7][1]*P[13][0]*P[15][0]*P[24][0] + 
   P[2][0]*P[6][0]*P[13][1]*P[15][0]*P[24][0] - 
   P[1][0]*P[7][0]*P[13][1]*P[15][0]*P[24][0] + 
   P[3][0]*P[7][0]*P[11][0]*P[15][1]*P[24][0] - 
   P[2][0]*P[8][0]*P[11][0]*P[15][1]*P[24][0] - 
   P[3][0]*P[6][0]*P[12][0]*P[15][1]*P[24][0] + 
   P[1][0]*P[8][0]*P[12][0]*P[15][1]*P[24][0] + 
   P[2][0]*P[6][0]*P[13][0]*P[15][1]*P[24][0] - 
   P[1][0]*P[7][0]*P[13][0]*P[15][1]*P[24][0] - 
   P[3][1]*P[7][0]*P[10][0]*P[16][0]*P[24][0] - 
   P[3][0]*P[7][1]*P[10][0]*P[16][0]*P[24][0] + 
   P[2][1]*P[8][0]*P[10][0]*P[16][0]*P[24][0] + 
   P[2][0]*P[8][1]*P[10][0]*P[16][0]*P[24][0] - 
   P[3][0]*P[7][0]*P[10][1]*P[16][0]*P[24][0] + 
   P[2][0]*P[8][0]*P[10][1]*P[16][0]*P[24][0] + 
   P[3][1]*P[5][0]*P[12][0]*P[16][0]*P[24][0] + 
   P[3][0]*P[5][1]*P[12][0]*P[16][0]*P[24][0] - 
   P[0][1]*P[8][0]*P[12][0]*P[16][0]*P[24][0] - 
   P[0][0]*P[8][1]*P[12][0]*P[16][0]*P[24][0] + 
   P[3][0]*P[5][0]*P[12][1]*P[16][0]*P[24][0] - 
   P[0][0]*P[8][0]*P[12][1]*P[16][0]*P[24][0] - 
   P[2][1]*P[5][0]*P[13][0]*P[16][0]*P[24][0] - 
   P[2][0]*P[5][1]*P[13][0]*P[16][0]*P[24][0] + 
   P[0][1]*P[7][0]*P[13][0]*P[16][0]*P[24][0] + 
   P[0][0]*P[7][1]*P[13][0]*P[16][0]*P[24][0] - 
   P[2][0]*P[5][0]*P[13][1]*P[16][0]*P[24][0] + 
   P[0][0]*P[7][0]*P[13][1]*P[16][0]*P[24][0] - 
   P[3][0]*P[7][0]*P[10][0]*P[16][1]*P[24][0] + 
   P[2][0]*P[8][0]*P[10][0]*P[16][1]*P[24][0] + 
   P[3][0]*P[5][0]*P[12][0]*P[16][1]*P[24][0] - 
   P[0][0]*P[8][0]*P[12][0]*P[16][1]*P[24][0] - 
   P[2][0]*P[5][0]*P[13][0]*P[16][1]*P[24][0] + 
   P[0][0]*P[7][0]*P[13][0]*P[16][1]*P[24][0] + 
   P[3][1]*P[6][0]*P[10][0]*P[17][0]*P[24][0] + 
   P[3][0]*P[6][1]*P[10][0]*P[17][0]*P[24][0] - 
   P[1][1]*P[8][0]*P[10][0]*P[17][0]*P[24][0] - 
   P[1][0]*P[8][1]*P[10][0]*P[17][0]*P[24][0] + 
   P[3][0]*P[6][0]*P[10][1]*P[17][0]*P[24][0] - 
   P[1][0]*P[8][0]*P[10][1]*P[17][0]*P[24][0] - 
   P[3][1]*P[5][0]*P[11][0]*P[17][0]*P[24][0] - 
   P[3][0]*P[5][1]*P[11][0]*P[17][0]*P[24][0] + 
   P[0][1]*P[8][0]*P[11][0]*P[17][0]*P[24][0] + 
   P[0][0]*P[8][1]*P[11][0]*P[17][0]*P[24][0] - 
   P[3][0]*P[5][0]*P[11][1]*P[17][0]*P[24][0] + 
   P[0][0]*P[8][0]*P[11][1]*P[17][0]*P[24][0] + 
   P[1][1]*P[5][0]*P[13][0]*P[17][0]*P[24][0] + 
   P[1][0]*P[5][1]*P[13][0]*P[17][0]*P[24][0] - 
   P[0][1]*P[6][0]*P[13][0]*P[17][0]*P[24][0] - 
   P[0][0]*P[6][1]*P[13][0]*P[17][0]*P[24][0] + 
   P[1][0]*P[5][0]*P[13][1]*P[17][0]*P[24][0] - 
   P[0][0]*P[6][0]*P[13][1]*P[17][0]*P[24][0] + 
   P[3][0]*P[6][0]*P[10][0]*P[17][1]*P[24][0] - 
   P[1][0]*P[8][0]*P[10][0]*P[17][1]*P[24][0] - 
   P[3][0]*P[5][0]*P[11][0]*P[17][1]*P[24][0] + 
   P[0][0]*P[8][0]*P[11][0]*P[17][1]*P[24][0] + 
   P[1][0]*P[5][0]*P[13][0]*P[17][1]*P[24][0] - 
   P[0][0]*P[6][0]*P[13][0]*P[17][1]*P[24][0] - 
   P[2][1]*P[6][0]*P[10][0]*P[18][0]*P[24][0] - 
   P[2][0]*P[6][1]*P[10][0]*P[18][0]*P[24][0] + 
   P[1][1]*P[7][0]*P[10][0]*P[18][0]*P[24][0] + 
   P[1][0]*P[7][1]*P[10][0]*P[18][0]*P[24][0] - 
   P[2][0]*P[6][0]*P[10][1]*P[18][0]*P[24][0] + 
   P[1][0]*P[7][0]*P[10][1]*P[18][0]*P[24][0] + 
   P[2][1]*P[5][0]*P[11][0]*P[18][0]*P[24][0] + 
   P[2][0]*P[5][1]*P[11][0]*P[18][0]*P[24][0] - 
   P[0][1]*P[7][0]*P[11][0]*P[18][0]*P[24][0] - 
   P[0][0]*P[7][1]*P[11][0]*P[18][0]*P[24][0] + 
   P[2][0]*P[5][0]*P[11][1]*P[18][0]*P[24][0] - 
   P[0][0]*P[7][0]*P[11][1]*P[18][0]*P[24][0] - 
   P[1][1]*P[5][0]*P[12][0]*P[18][0]*P[24][0] - 
   P[1][0]*P[5][1]*P[12][0]*P[18][0]*P[24][0] + 
   P[0][1]*P[6][0]*P[12][0]*P[18][0]*P[24][0] + 
   P[0][0]*P[6][1]*P[12][0]*P[18][0]*P[24][0] - 
   P[1][0]*P[5][0]*P[12][1]*P[18][0]*P[24][0] + 
   P[0][0]*P[6][0]*P[12][1]*P[18][0]*P[24][0] - 
   P[2][0]*P[6][0]*P[10][0]*P[18][1]*P[24][0] + 
   P[1][0]*P[7][0]*P[10][0]*P[18][1]*P[24][0] + 
   P[2][0]*P[5][0]*P[11][0]*P[18][1]*P[24][0] - 
   P[0][0]*P[7][0]*P[11][0]*P[18][1]*P[24][0] - 
   P[1][0]*P[5][0]*P[12][0]*P[18][1]*P[24][0] + 
   P[0][0]*P[6][0]*P[12][0]*P[18][1]*P[24][0] + 
   P[3][0]*P[7][0]*P[11][0]*P[15][0]*P[24][1] - 
   P[2][0]*P[8][0]*P[11][0]*P[15][0]*P[24][1] - 
   P[3][0]*P[6][0]*P[12][0]*P[15][0]*P[24][1] + 
   P[1][0]*P[8][0]*P[12][0]*P[15][0]*P[24][1] + 
   P[2][0]*P[6][0]*P[13][0]*P[15][0]*P[24][1] - 
   P[1][0]*P[7][0]*P[13][0]*P[15][0]*P[24][1] - 
   P[3][0]*P[7][0]*P[10][0]*P[16][0]*P[24][1] + 
   P[2][0]*P[8][0]*P[10][0]*P[16][0]*P[24][1] + 
   P[3][0]*P[5][0]*P[12][0]*P[16][0]*P[24][1] - 
   P[0][0]*P[8][0]*P[12][0]*P[16][0]*P[24][1] - 
   P[2][0]*P[5][0]*P[13][0]*P[16][0]*P[24][1] + 
   P[0][0]*P[7][0]*P[13][0]*P[16][0]*P[24][1] + 
   P[3][0]*P[6][0]*P[10][0]*P[17][0]*P[24][1] - 
   P[1][0]*P[8][0]*P[10][0]*P[17][0]*P[24][1] - 
   P[3][0]*P[5][0]*P[11][0]*P[17][0]*P[24][1] + 
   P[0][0]*P[8][0]*P[11][0]*P[17][0]*P[24][1] + 
   P[1][0]*P[5][0]*P[13][0]*P[17][0]*P[24][1] - 
   P[0][0]*P[6][0]*P[13][0]*P[17][0]*P[24][1] - 
   P[2][0]*P[6][0]*P[10][0]*P[18][0]*P[24][1] + 
   P[1][0]*P[7][0]*P[10][0]*P[18][0]*P[24][1] + 
   P[2][0]*P[5][0]*P[11][0]*P[18][0]*P[24][1] - 
   P[0][0]*P[7][0]*P[11][0]*P[18][0]*P[24][1] - 
   P[1][0]*P[5][0]*P[12][0]*P[18][0]*P[24][1] + 
				      P[0][0]*P[6][0]*P[12][0]*P[18][0]*P[24][1]
};
  };

#endif
  
}//namespace