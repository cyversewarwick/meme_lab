wget http://search.cpan.org/CPAN/authors/id/L/LD/LDS/GD-2.56.tar.gz
tar xvzf GD-2.56.tar.gz
cd GD*
perl -i~ -pE'say "Getopt::Long::Configure(qw( pass_through ));" if /GetOptions/;' Build.PL
/usr/bin/perl Build.PL --installdirs site
chmod 755 Build.PL
./Build.PL installdeps
./Build installdeps
./Build test
./Build install