function [urlConnection,errorid,errormsg] = urlreadwrite(fcn,urlChar)
%URLREADWRITE A helper function for URLREAD and URLWRITE.

%   Matthew J. Simoneau, June 2005
%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.3.6.1 $ $Date: 2009/01/30 22:37:42 $


%Copyright (c) 2010, Dan Ellis
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are
%met:
%
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in
%      the documentation and/or other materials provided with the distribution
%    * Neither the name of the Columbia University nor the names
%      of its contributors may be used to endorse or promote products derived
%      from this software without specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%POSSIBILITY OF SUCH DAMAGE.



% Default output arguments.
urlConnection = [];
errorid = '';
errormsg = '';

% Determine the protocol (before the ":").
protocol = urlChar(1:min(find(urlChar==':'))-1);

% Try to use the native handler, not the ice.* classes.
switch protocol
    case 'http'
        try
            handler = sun.net.www.protocol.http.Handler;
        catch exception %#ok
            handler = [];
        end
    case 'https'
        try
            handler = sun.net.www.protocol.https.Handler;
        catch exception %#ok
            handler = [];
        end
    otherwise
        handler = [];
end

% Create the URL object.
try
    if isempty(handler)
        url = java.net.URL(urlChar);
    else
        url = java.net.URL([],urlChar,handler);
    end
catch exception %#ok
    errorid = ['MATLAB:' fcn ':InvalidUrl'];
    errormsg = 'Either this URL could not be parsed or the protocol is not supported.';
    return
end

% Get the proxy information using MathWorks facilities for unified proxy
% prefence settings.
mwtcp = com.mathworks.net.transport.MWTransportClientPropertiesFactory.create();
proxy = mwtcp.getProxy(); 


% Open a connection to the URL.
if isempty(proxy)
    urlConnection = url.openConnection;
else
    urlConnection = url.openConnection(proxy);
end
