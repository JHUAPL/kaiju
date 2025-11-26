Setting up and sharing NCAR GLADE data via a Globus Collections
===============================================================

Users with accounts on NCAR's HPC systems can now share data with other users
via Globus.  The other users will not need an NCAR account, but depending on
the permissions you select for sharing they may need to setup a free Globus
account to access the data.

First a bit of terminology.  Users interact with collections and a *mapped
collection* is a set of files hosted at an endpoint. The NCAR GLADE filesystem
is a mapped collection. Users can create a guest collection from a mapped
collection with distinct sharing permissions.

You will need to authenticate to the *NCAR mapped collections* via
`Globus.org <https://app.globus.org/>`_ with your Globus Id.  Using the File
manager option you will need search for NCAR GLADE collection.  Choose the one
with the description for GridFTP access to GLADE filesystem using UCAS token
authentication.  You will need to login with your NCAR id and duo token. If
you want to share data from the NCAR campaign storage system you need to use
the file manager to search for NCAR Campaign Storage.  Remember that the
campaign storage system is not accessible from cheyenne and data must be
transfer to it via casper.

It is highly recommended to create bookmarks of the NCAR Glade and NCAR
Campaign collections as they will be the starting point for creating any guest
collection to share data with other users.

.. image:: https://bitbucket.org/repo/kMoBzBp/images/551980672-GlobusBookmarks.png
   :target: https://bitbucket.org/repo/kMoBzBp/images/551980672-GlobusBookmarks.png
   :alt: GlobusBookmarks.png


To share a directory with other Globus users first create the directory on the
NCAR filesystem you wish to share with other users.  Via the Globus website
use your bookmarks to navigate to the file system you just created and then
select the share option to create a new guest collection.  The image below
show how to share the directory
/glade/campaign/hao/msphere/wiltbemj/ExampleShare as a guest collection.


.. image:: https://bitbucket.org/repo/kMoBzBp/images/2261017910-GlobusShare.png
   :target: https://bitbucket.org/repo/kMoBzBp/images/2261017910-GlobusShare.png
   :alt: GlobusShare.png


Once you click share you will need to choose the option for creating a new
guest collection and then you'll need to provide the required display name and
meta data.  If the data is related to a publication you can provide the DOI in
the information link available in *view more fields* option on the web page.

Once you create the collection you will have option to add permissions for
sharing the data with other users.    There are currently for levels of
sharing, specific users, groups, all global users, and public.  The web site
will then provide a
`link  <https://app.globus.org/file-manager?origin_id=92ac5357-f5e6-4b01-bbf4-b1bb8c06af1f&origin_path=%2F>`_
that you can use for sharing with other users.
