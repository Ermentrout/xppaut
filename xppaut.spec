Name:           xppaut
Version:        8.0
Release:        1%{?dist}
Summary:        X-windows Phase Plane Plus Auto

Group:          Applications/Engineering
License:        GPL
URL:            http://www.math.pitt.edu/~bard/xpp/xpp.html
BuildRoot:      %{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)
Source0:        %{name}-%{version}.tar.gz
BuildRequires:  gcc libX11-devel

%define debug_package %{nil}

%description
XPPAUT is a tool for solving differential equations, difference equations, delay equations,
functional equations, boundary value problems, and stochastic equations.

%prep
%setup -q
%build
make %{?_smp_mflags}

%install
rm -rf %{buildroot}
rm -f ode/*.so
make BINDIR=%{_bindir} MANDIR=%{_mandir}/man1 DESTDIR=%{buildroot} install
%clean
rm -rf %{buildroot}
%files
%defattr(-,root,root,-)
/usr/bin/xppaut
%doc %attr(0444,root,root) /usr/share/man/man1/xppaut.1.gz
%doc %attr(0444,root,root) /usr/share/doc/xppaut/*
%changelog
