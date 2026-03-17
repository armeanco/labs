std::string alphabetWar(const std::string& fight) {
  int cnt_l = 0, cnt_r = 0;
  for(int i = 0; i < static_cast<int>(fight.size()); ++i) {
    if(fight[i] == 'w' && fight[i - 1] != '*' && fight[i + 1] != '*') cnt_l += 4;
    if(fight[i] == 'p' && fight[i - 1] != '*' && fight[i + 1] != '*') cnt_l += 3;
    if(fight[i] == 'b' && fight[i - 1] != '*' && fight[i + 1] != '*') cnt_l += 2;
    if(fight[i] == 's' && fight[i - 1] != '*' && fight[i + 1] != '*') cnt_l += 1;
    if(fight[i] == 'm' && fight[i - 1] != '*' && fight[i + 1] != '*') cnt_r += 4;
    if(fight[i] == 'q' && fight[i - 1] != '*' && fight[i + 1] != '*') cnt_r += 3;
    if(fight[i] == 'd' && fight[i - 1] != '*' && fight[i + 1] != '*') cnt_r += 2;
    if(fight[i] == 'z' && fight[i - 1] != '*' && fight[i + 1] != '*') cnt_r += 1;
  }
  return cnt_l > cnt_r ? "Left side wins!" : cnt_l == cnt_r ? "Let's fight again!" : "Right side wins!";
}
