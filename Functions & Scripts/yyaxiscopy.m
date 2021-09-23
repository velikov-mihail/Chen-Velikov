function yyaxiscopy
% copy left axis onto right
% 2018 06 Andrew

yyaxis left
ylimtemp = ylim;
yticktemp = get(gca,'YTick');

yyaxis right
set(gca,'YColor','k')
ylim(ylimtemp)
set(gca,'YTick',yticktemp)

yyaxis left
