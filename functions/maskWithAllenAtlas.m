function [hits] = maskWithAllenAtlas(params,subs,maskids,skip)
% returns the list of inds that hits the Allen mask
addpath(genpath('/groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/matlab_dbqueries'))
addpath(genpath('/groups/mousebrainmicro/mousebrainmicro/Software/Matlab/Registration tools'))
allenH5 = '/groups/mousebrainmicro/mousebrainmicro/Software/Matlab/Registration tools/Dependencies/OntologyAtlas.h5';
transformFile = '/groups/mousebrainmicro/mousebrainmicro/registration/Database/2018-08-01/Transform.2018-08-01.h5';

if nargin==1
    maskids = [0 73 81 89 98 108 116 124 129 140 145 153 164]; % outside and ventrical system
    maskids = [56 672 1022 1031]; %ACB/Caudoputamen/Globus pallidus, external segment/Globus pallidus, internal segment
    maskids = [956 844 882 686 672 56 1022 1031 1021 1085 719 882 583 182305705 182305709 182305713]
end
%%
% out = transSample2Allen(sampleArg, transformFile,outputFolder)
try
    tic
    vectorField = h5read(transformFile,'/DisplacementField');
    toc
    tFormVecField = h5readatt(transformFile,'/DisplacementField','Transformation_Matrix');
catch
    error('Could not read Displacement Field file: %s',transformFile);
end

%% voxel -> um -> allen
A_vox2um = [[diag([params.voxres]) [params.ox params.oy params.oz]'/1e3];[0 0 0 1]];
subs_um = [subs ones(size(subs,1),1)]*A_vox2um';
subs_allen = subs_um*tFormVecField';

%%
pixPosVecField = ceil(subs_allen);
nPoints = size(pixPosVecField,1);

% ravel multiindex
dims = size(vectorField);
inds = sub2ind(dims(1:3),pixPosVecField(:,1),pixPosVecField(:,2),pixPosVecField(:,3));
rvectorField = reshape(vectorField,[],3);
vector = rvectorField(inds,:);
subs_um_updated = ceil(subs_um(:,1:3) + vector);

%%
allentForm = h5readatt(allenH5,'/OntologyAtlas','Transformation_Matrix');
atlas = h5read(allenH5,'/OntologyAtlas');

mask = 0*atlas;
for ii = 1:length(maskids)
    mask=mask|atlas==maskids(ii);
end
mask = imdilate(mask,strel('sphere',10));
%%
dims_atlas = size(atlas);
pixPos = ceil([subs_um_updated,ones(nPoints,1)]*allentForm');
val_outofbounds = (pixPos(:,1)>dims_atlas(:,1) | pixPos(:,2)>dims_atlas(:,2) | pixPos(:,3)>dims_atlas(:,3) | any(pixPos<1,2));
inds_inofbounds = find(~val_outofbounds);
pixPos_ = pixPos(inds_inofbounds,:);
indPoints = sub2ind(dims_atlas,pixPos_(:,1),pixPos_(:,2),pixPos_(:,3));
hits = inds_inofbounds(mask(indPoints));

return
%%
clc
for ii=1:length(fg.allenMesh);
    if strfind(fg.allenMesh(ii).hierarchy_path,'/507/')
        ii
        fg.allenMesh(ii).hierarchy_path
        fg.allenMesh(ii)
%         break;
    end;
    ii;
end


%%
junk = indPoints(mask(indPoints));
theseinds = find(~outofbounds);
subsjunk=[];
[subsjunk(:,1),subsjunk(:,2),subsjunk(:,3)] = ind2sub(dims_atlas,junk);

subs_ = subs;
subs_(junk,:)=[];
% bank = nan(nPoints,1);
% bank(theseinds) = atlas(indPoints);
% bank(bank==0)=nan;
%%
compartmentlist = [997 8 567 688 695 315 184 68 667 500 107 219 299 644 947 985 320 943 648 844 882 993 656 962 767 1021 1085 453 12993 12994 12995 12996 12997 12998 322 793 346 865 921 686 719 353 558 838 654 702 889 929 329 981 201 1047 1070 1038 1062 337 1030 113 1094 1128 478 510 345 878 657 950 974 1102 2 369 450 854 577 625 945 1026 361 1006 670 1086 1111 9 461 182305689 182305693 182305697 182305701 182305705 182305709 182305713 378 873 806 1035 1090 862 893 1057 36 180 148 187 638 662 677 897 1106 1010 1058 857 849 247 1011 527 600 678 252 156 243 1002 735 251 816 847 954 1005 1027 696 643 759 791 249 456 1018 959 755 990 1023 520 598 669 801 561 913 937 457 497 402 1074 905 1114 233 601 649 394 281 1066 401 433 1046 441 409 421 973 573 613 74 121 385 593 821 721 778 33 305 425 750 269 869 902 377 393 533 805 41 501 565 257 469 31 572 1053 739 179 227 39 935 211 1015 919 927 48 588 296 772 810 819 972 171 195 304 363 84 132 44 707 747 556 827 1054 1081 714 264 492 352 476 516 723 448 412 630 440 488 731 484 524 582 620 910 738 746 969 288 1125 608 680 95 104 996 328 1101 783 831 111 120 163 344 314 355 119 704 694 800 675 699 254 894 671 965 774 906 279 879 442 434 545 610 274 330 886 542 606 430 687 590 622 22 532 241 635 683 308 340 541 97 1127 234 289 729 786 922 335 368 540 692 888 895 836 427 988 977 1045 698 507 212 220 228 236 244 151 188 196 204 159 167 175 183 191 199 160 168 589 597 297 1034 1042 1050 1059 605 306 1067 1075 1082 814 496 535 360 646 267 961 152 276 284 291 619 392 260 268 1139 631 639 192 200 208 647 655 584 376 216 224 232 663 592 383 240 248 256 788 400 408 416 424 566 517 1140 1141 1142 1089 1080 375 382 391 399 407 415 423 431 438 446 454 463 471 479 486 495 504 726 10703 10704 632 10702 734 742 751 758 766 775 782 790 799 807 815 823 982 19 822 909 918 1121 20 999 715 764 52 92 312 139 387 28 60 926 526 543 468 508 664 712 727 550 743 934 259 324 371 419 1133 843 10693 10694 10695 1037 10696 10697 10698 1084 10699 10700 10701 502 509 829 845 837 518 853 870 861 703 16 583 942 952 966 131 295 303 311 451 319 327 334 780 623 477 485 672 493 56 998 754 481 489 144 458 465 473 275 242 250 258 266 310 333 278 23 292 536 544 551 559 1105 403 411 418 426 472 480 487 435 803 818 1022 1031 835 342 298 826 904 564 596 581 809 351 359 537 498 505 513 546 521 554 562 529 367 569 578 585 594 602 287 343 1129 549 864 637 629 685 709 718 725 733 741 406 414 422 609 1044 1008 475 1072 1079 1088 170 856 138 218 1020 1029 325 239 255 127 1096 1104 64 1120 1113 155 444 59 362 617 626 636 366 1077 571 149 15 181 51 189 599 907 575 930 262 1014 27 178 300 316 321 958 483 186 953 1097 157 390 332 432 38 71 47 79 103 652 660 94 55 87 110 30 118 223 141 72 80 263 272 830 668 676 684 452 523 763 914 1109 1124 126 133 347 286 338 689 467 88 700 708 716 724 331 210 491 732 525 1110 1118 557 1126 1 515 740 748 756 980 1004 63 439 447 455 464 693 761 769 777 785 946 290 194 226 356 364 173 470 614 797 796 804 10671 313 339 302 851 842 834 4 811 820 828 580 271 874 460 323 381 749 246 128 539 548 555 294 26 42 17 10 494 503 511 795 50 67 587 1100 215 531 628 634 706 1061 616 214 35 975 115 757 231 66 75 58 615 348 374 1052 165 12 100 197 591 872 1065 771 1132 612 82 90 99 7 867 123 881 860 868 875 883 891 890 899 915 923 398 122 105 114 987 280 880 283 898 931 1093 552 318 462 534 574 621 1117 679 137 130 147 162 604 146 238 350 358 354 386 207 607 112 560 96 101 720 711 1039 903 642 651 659 666 674 682 691 429 437 445 77 53 61 45 69 789 370 653 568 661 576 640 135 939 143 839 887 1048 372 83 136 106 203 235 955 963 307 395 1098 1107 852 859 938 970 978 154 161 177 169 995 1069 185 193 701 209 202 225 217 765 773 781 76 379 206 230 222 512 528 645 912 10707 10706 10705 920 976 10710 10709 10708 984 10713 10712 10711 928 992 10716 10715 10714 1001 10719 10718 10717 1091 10722 10721 10720 936 10725 10724 10723 944 10728 10727 10726 951 10731 10730 10729 957 10734 10733 10732 968 10737 10736 10735 1073 1007 10674 10673 10672 1017 1056 10677 10676 10675 1064 10680 10679 10678 1025 10683 10682 10681 1033 10686 10685 10684 1041 10689 10688 10687 1049 10692 10691 10690 1144 1145 1143 519 989 91 846 1009 967 885 949 840 1016 21 665 538 459 900 848 876 916 336 117 125 357 832 62 158 911 384 710 901 93 229 705 794 798 1131 1116 933 1076 413 948 841 641 506 658 633 482 808 917 237 717 813 925 792 932 570 522 858 586 514 380 388 396 697 871 29 389 245 261 270 293 277 253 285 627 960 744 752 326 812 85 850 866 78 1123 553 499 650 490 404 410 373 728 983 776 956 579 964 1108 971 979 986 784 6 924 1036 1012 1003 994 190 198 1019 1028 896 1092 14 86 365 1000 760 142 102 109 134 309 317 877 1051 1060 1043 863 397 221 736 855 205 213 941 991 768 884 892 908 940 1099 466 530 603 745 420 737 428 436 618 443 449 713 474 37 301 824 54 405 174 349 817 825 833 166 341 182 762 770 779 787 150 46 753 690 681 673 1068 722 1083 802 595 611 730 70 547 563 73 81 89 98 108 116 124 129 140 145 153 164 1024 1032 1055 1063 1071 1078 1040 1087 1095 1103 1112 1119 3 11 18 25 34 43 49 57 65 624 304325711];

%%
mesh = load('/nrs/mouselight/Shared Files/Mesh Info/allenMesh.mat');
ventMesh = mesh.allenMesh(find([mesh.allenMesh.id]==73));
in = intriangulation(ventMesh.v,ventMesh.f,pixPos(testThese,1:3));

%%
out = [];
for ii=1:length(compartmentlist)
    ii=10674
    if ~isnan(bank(iPnt))
        out{ii} = webread(sprintf('http://api.brain-map.org/api/v2/data/Structure/query.json?criteria=[id$eq%i]'...
            ,ii));
    else
        %                 fprintf('\n%f, %f, %f\n\tID: NaN, Outside',bank(iPnt,1));
    end
end
%%

[out,varargout] = allenCoord2allenIdSilent(subs_um_updated);



% %%
% vector = zeros(nPoints,3);
% for iPoint = 1:nPoints
%     vector(iPoint,:) = reshape(vectorField(pixPosVecField(iPoint,1),pixPosVecField(iPoint,2),...
%         pixPosVecField(iPoint,3),:),1,3);
% end


%%
% Apply transform.
fprintf('\nTransforming coordinates');
for iCell = 1:size(swcData,2)
    nPoints = size(swcData(iCell).Points,1);
    pixPosVecField = ceil(tFormVecField*[swcData(iCell).Points,ones(nPoints,1)]')';
    vector = zeros(nPoints,3);
    try
        for iPoint = 1:nPoints
            vector(iPoint,:) = reshape(vectorField(pixPosVecField(iPoint,1),pixPosVecField(iPoint,2),...
                pixPosVecField(iPoint,3),:),1,3);
        end
    catch
        error('Location(s) out-of-bounds for: %s',swcData(iCell).Name);
    end
    % register.
    swcData(iCell).Points = swcData(iCell).Points + vector;
    swcData(iCell).Swc(:,3:5) = swcData(iCell).Points;
end