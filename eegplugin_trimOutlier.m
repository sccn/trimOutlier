% eegplugin_trimOutlier() - calls pop_trimOutlier

% History:
% 03/07/2014 ver 1.0 by Makoto. Former firstpassOutlierTrimmer redesigned and renamed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 05/16/2013 ver 1.0 by Makoto. Created.

% Author: Makoto Miyakoshi, SCCN,INC,UCSD 2013
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function eegplugin_trimOutlier( fig, try_strings, catch_strings)

vers = '1.0';
if nargin < 3
    error('eegplugin_trimOutlier requires 3 arguments');
end

% create menu
toolsmenu = findobj(fig, 'tag', 'tools');
uimenu( toolsmenu, 'label', 'trimOutlier', 'separator','on',...
    'callback', 'EEG = pop_trimOutlier(EEG); [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); eeglab redraw');
