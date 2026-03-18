int comparePowers(std::pair<long, long> n1, std::pair<long, long> n2){
  int _ = 0, $ = 0;
  while(n1.first >>= 1) _++;
  while(n2.first >>= 1) $++;
  return n1.second *_> n2.second *$ ? -1 : n1.second *_< n2.second *$ ? 1 : 0;
}
