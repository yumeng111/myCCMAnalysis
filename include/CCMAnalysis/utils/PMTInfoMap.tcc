#ifndef PMTInfoMap_tcc
#define PMTInfoMap_tcc

/*!************************************************************************************************
 * \fn void PMTInfoMap::DecodeKey(T key, T & digit, T & channel)
 * \brief Static function that decodes the map key
 * \param[in] key The key value
 * \param[out] digit Board number
 * \param[out] channel Channel number
 **************************************************************************************************/
template<typename T, typename U>
void PMTInfoMap::DecodeKey(T key, U & digit, U & channel) {
  digit = key/16;
  channel = key%16;

  return;
}

#endif // PMTInfoMap_h
